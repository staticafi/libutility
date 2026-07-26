#ifndef VISUALIZER_BREAKPOINT_HPP_INCLUDED
#   define VISUALIZER_BREAKPOINT_HPP_INCLUDED

#   include <cstdint>

#   define VISUALIZER_BREAKPOINT()                                                          \
        do { static ::visualizer::detail::BreakPointID const id =                           \
                ::visualizer::detail::register_breakpoint(__FILE__,__LINE__,__FUNCTION__);  \
             ::visualizer::detail::visualize(id); } while (false)

namespace visualizer::detail {

    struct BreakPoint
    {
        char const* file{ nullptr };
        int line{ 0 };
        char const* func{ nullptr };
    };

    using BreakPointID = std::uint32_t;

    BreakPointID register_breakpoint(char const* file, int line, char const* func);
    BreakPoint const* get_breakpoint(BreakPointID id);
    BreakPointID get_current_breakpoint_id();
    BreakPointID get_end_breakpoint_id(); // The count of all breakpoints. Valid IDs are in range [0, get_end_breakpoint_id()).
    void  visualize(BreakPointID id);
}

#endif
