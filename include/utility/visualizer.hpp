#ifndef VISUALIZER_HPP_INCLUDED
#   define VISUALIZER_HPP_INCLUDED

#   include <functional>
#   include <memory>
#   include <mutex>

namespace visualizer {


struct VisualizerBase
{
    static std::mutex  s_mutex;
    static bool s_stop_flag;
    static bool s_render;

    VisualizerBase() {}
    virtual ~VisualizerBase() {}

    virtual void next_frame() {}

    bool can_render() const;
    void set_waiting_for_content() const;

private:

    VisualizerBase(VisualizerBase const&) = delete;
    VisualizerBase(VisualizerBase&&) = delete;
    VisualizerBase& operator=(VisualizerBase const&) = delete;
    VisualizerBase& operator=(VisualizerBase&&) = delete;
};

using ConstructorType = std::function<std::unique_ptr<VisualizerBase>()>;

void  create_visualizer(ConstructorType const& constructor);
void  destroy_visualizer();

struct  VisualizerTerminationException {};


}

#endif
