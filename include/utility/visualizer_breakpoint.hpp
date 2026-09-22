#ifndef VISUALIZER_BREAKPOINT_HPP_INCLUDED
#   define VISUALIZER_BREAKPOINT_HPP_INCLUDED

#   include <source_location>
#   include <cstdint>

#   define VISUALIZER_BREAKPOINT()                                                          \
        do { static ::visualizer::BreakPointID const id =                                   \
                ::visualizer::register_breakpoint(                                          \
                    std::source_location::current().file_name(),                            \
                    std::source_location::current().line(),                                 \
                    std::source_location::current().function_name()                         \
                    );                                                                      \
             ::visualizer::on_breakpoint_hit(id); } while (false)

namespace visualizer {

    struct BreakPoint
    {
        char const* file{ nullptr };
        int line{ 0 };
        char const* func{ nullptr };
    };

    // The ID 0 is special. It represent the nearest breakpoint
    // to be hit by program execution from the current breakpoint.
    // The other IDs (greater than 0) are concrete breakpoints.
    using BreakPointID = std::uint32_t;

    // Returns ID greater than zero, because ID 0 is reserved for the special breakpoint.
    BreakPointID register_breakpoint(char const* file, int line, char const* func);

    // The count excludes the special breakpoint with ID 0.
    std::uint32_t get_num_registered_breakpoints();

    // The passed ID must be greater than zero.
    BreakPoint const* get_registered_breakpoint(BreakPointID id);

    BreakPointID get_target_breakpoint_id();
    void set_target_breakpoint_id(BreakPointID id);

    bool is_execution_paused_on_breakpoint();
    void request_resume_execution();

    void  on_breakpoint_hit(BreakPointID id);

}

#endif
