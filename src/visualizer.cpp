#include <utility/visualizer.hpp>
#include <utility/visualizer_breakpoint.hpp>
#include <thread>
#include <chrono>
#include <memory>
#include <stdexcept>
#include <iostream>

namespace visualizer {


std::mutex  VisualizerBase::s_mutex{};
bool VisualizerBase::s_stop_flag{ false };
bool VisualizerBase::s_render{ false };


bool VisualizerBase::can_render() const
{
    bool do_render;
    {
        std::lock_guard<std::mutex> const lock(s_mutex);
        do_render = !s_stop_flag && s_render;
    }
    return do_render;
}


void VisualizerBase::set_waiting_for_content() const
{
    std::lock_guard<std::mutex> const lock(s_mutex);
    s_render = false;
}


///////////////////////////////////////////////////////////////////////////
// Next follows thread managements and interaction between threads.
///////////////////////////////////////////////////////////////////////////


static std::thread  s_visualizer_thread{};


static void visualizer_thread_procedure(ConstructorType const& constructor)
{
    std::unique_ptr<VisualizerBase>  visualizer_ptr{ constructor() };
    while (true)
    {
        {
            std::lock_guard<std::mutex> const lock(VisualizerBase::s_mutex);
            if (VisualizerBase::s_stop_flag)
                break;
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
    std::lock_guard<std::mutex> const lock(VisualizerBase::s_mutex);
    visualizer_ptr = nullptr;
    VisualizerBase::s_stop_flag = true;
}


void  create_visualizer(ConstructorType const& constructor)
{
    if (s_visualizer_thread.joinable())
        return;
    VisualizerBase::s_stop_flag = false;
    VisualizerBase::s_render = false;
    s_visualizer_thread = std::thread(visualizer_thread_procedure, constructor);
}


void  destroy_visualizer()
{
    if (!s_visualizer_thread.joinable())
        return;

    {
        std::lock_guard<std::mutex> const lock(VisualizerBase::s_mutex);
        VisualizerBase::s_stop_flag = true;
        VisualizerBase::s_render = false;
    }
    if (s_visualizer_thread.joinable())
        s_visualizer_thread.join();
}


}

namespace visualizer::detail {


static std::vector<BreakPoint> s_breakpoints{};
static BreakPointID s_current_breakpoint_id{ 0U };


BreakPointID register_breakpoint(char const* const file, int const line, char const* const func)
{
    s_breakpoints.push_back(BreakPoint{
        .file = file,
        .line = line,
        .func = func
    });
    return (BreakPointID)(s_breakpoints.size() - 1ULL);
}


BreakPoint const*  get_breakpoint(BreakPointID const id)
{
    return &s_breakpoints.at(id);
}


BreakPointID get_current_breakpoint_id()
{
    return s_current_breakpoint_id;
}


BreakPointID get_end_breakpoint_id()
{
    return (BreakPointID)s_breakpoints.size();
}


void  visualize(BreakPointID const id)
{
    if (!s_visualizer_thread.joinable())
        return;

    {
        std::lock_guard<std::mutex> const lock(VisualizerBase::s_mutex);
        if (VisualizerBase::s_stop_flag)
            throw VisualizerTerminationException{};
        s_current_breakpoint_id = id;
        VisualizerBase::s_render = true;
    }

    while (s_visualizer_thread.joinable())
    {
        {
            std::lock_guard<std::mutex> const lock(VisualizerBase::s_mutex);
            if (VisualizerBase::s_stop_flag)
                throw VisualizerTerminationException{};
            if (!VisualizerBase::s_render)
                break;
        }
        std::this_thread::yield();
        std::this_thread::sleep_for(std::chrono::microseconds(100));
    }
}


}
