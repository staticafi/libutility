#ifndef VISUALIZER_HPP_INCLUDED
#   define VISUALIZER_HPP_INCLUDED

#   include <functional>
#   include <memory>

namespace visualizer {


struct Visualizer;
using ConstructorType = std::function<std::unique_ptr<Visualizer>()>;


struct Visualizer
{
    Visualizer() {}
    virtual ~Visualizer() {}

    virtual void next_frame() {}
    virtual void on_data_changed() {}

    void stop();

private:

    Visualizer(Visualizer const&) = delete;
    Visualizer(Visualizer&&) = delete;
    Visualizer& operator=(Visualizer const&) = delete;
    Visualizer& operator=(Visualizer&&) = delete;
};


struct  TerminationException {};


void  create(ConstructorType const& constructor);
void  destroy();


}

#endif
