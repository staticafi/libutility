#include <utility/visualizer.hpp>
#include <utility/visualizer_breakpoint.hpp>
#include <utility/assumptions.hpp>
#include <thread>
#include <chrono>
#include <mutex>
#include <stdexcept>
#include <iostream>

namespace visualizer {


enum struct ExecutionState
{
    BREAKPOINT_JUST_HIT,
    WAITING_ON_BREAKPOINT,
    RESUME_REQUESTED,
    RESUMED
};


static std::mutex  s_mutex{};
static std::thread  s_visualizer_thread{};
static bool  s_instance_created{ false };
static bool  s_stop_flag{ false };
static ExecutionState  s_execution_state{ ExecutionState::RESUMED };
static std::vector<BreakPoint> s_registered_breakpoints{};
static BreakPointID s_target_breakpoint{ 0U };


static void visualizer_thread_procedure(ConstructorType const& constructor)
{
    std::unique_ptr<Visualizer>  visualizer_ptr{ constructor() };
    std::this_thread::yield();
    std::this_thread::sleep_for(std::chrono::microseconds(100));
    while (true)
    {
        {
            std::lock_guard<std::mutex> const lock(s_mutex);
            if (s_stop_flag)
                break;
            if (s_execution_state == ExecutionState::BREAKPOINT_JUST_HIT)
            {
                visualizer_ptr->on_data_changed();
                s_execution_state = ExecutionState::WAITING_ON_BREAKPOINT;
            }
            else if (s_execution_state == ExecutionState::RESUME_REQUESTED)
                s_execution_state = ExecutionState::RESUMED;
        }
        try
        {
            visualizer_ptr->next_frame();
            //std::this_thread::yield();
            using namespace std::chrono_literals;
            std::this_thread::sleep_for(10ms);
        }
        catch (...)
        {
            break;
        }
    }
    std::lock_guard<std::mutex> const lock(s_mutex);
    visualizer_ptr = nullptr;
    s_stop_flag = true;
}


void  create(ConstructorType const& constructor)
{
    ASSUMPTION(!s_instance_created);
    if (s_visualizer_thread.joinable())
        return;
    s_stop_flag = false;
    s_execution_state = ExecutionState::RESUMED;
    s_visualizer_thread = std::thread(visualizer_thread_procedure, constructor);
    s_instance_created = true;
}


void  destroy()
{
    if (!s_visualizer_thread.joinable())
        return;

    {
        std::lock_guard<std::mutex> const lock(s_mutex);
        s_stop_flag = true;
        s_execution_state = ExecutionState::RESUMED;
    }

    if (s_visualizer_thread.joinable())
        s_visualizer_thread.join();
}


void Visualizer::stop()
{
    std::lock_guard<std::mutex> const lock(s_mutex);
    s_stop_flag = true;
    s_execution_state = ExecutionState::RESUMED;

}


BreakPointID register_breakpoint(char const* const file, int const line, char const* const func)
{
    s_registered_breakpoints.push_back(BreakPoint{
        .file = file,
        .line = line,
        .func = func
    });
    return (BreakPointID)s_registered_breakpoints.size();
}

std::uint32_t get_num_registered_breakpoints() { return (std::uint32_t)s_registered_breakpoints.size(); }

BreakPoint const* get_registered_breakpoint(BreakPointID const id)
{
    ASSUMPTION(id > 0U);
    return &s_registered_breakpoints.at(id - 1U);
}

BreakPointID get_target_breakpoint_id() { return s_target_breakpoint; }
void set_target_breakpoint_id(BreakPointID const id) { s_target_breakpoint = id; }


static bool  raw_execution_paused_condition()
{
    return s_instance_created && !s_stop_flag && s_execution_state == ExecutionState::WAITING_ON_BREAKPOINT;
}


bool is_execution_paused_on_breakpoint()
{
    bool result;
    {
        std::lock_guard<std::mutex> const lock(s_mutex);
        result = raw_execution_paused_condition();
    }
    return result;
}


void request_resume_execution()
{
    std::lock_guard<std::mutex> const lock(s_mutex);
    if (raw_execution_paused_condition())
        s_execution_state = ExecutionState::RESUME_REQUESTED;
}


void  on_breakpoint_hit(BreakPointID const id)
{
    if (!s_visualizer_thread.joinable())
        return;

    {
        std::lock_guard<std::mutex> const lock(s_mutex);
        if (s_stop_flag)
            throw TerminationException{};
        if (get_target_breakpoint_id() != 0U && id != get_target_breakpoint_id())
            return;
        set_target_breakpoint_id(id);
        s_execution_state = ExecutionState::BREAKPOINT_JUST_HIT;
    }

    while (s_visualizer_thread.joinable())
    {
        {
            std::lock_guard<std::mutex> const lock(s_mutex);
            if (s_stop_flag)
                throw TerminationException{};
            if (s_execution_state == ExecutionState::RESUMED)
                break;
        }
        std::this_thread::yield();
        std::this_thread::sleep_for(std::chrono::microseconds(100));
    }
}


}
