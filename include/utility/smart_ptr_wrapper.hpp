#ifndef UTILITY_SMART_PTR_WRAPPER_HPP_INCLUDED
#   define UTILITY_SMART_PTR_WRAPPER_HPP_INCLUDED

#   include <utility/config.hpp>
#   include <memory>

/// The only purpose of this module is to bypass the issue
/// with visualization of content of an STL container pointed
/// to by an std::shared/weak/unique pointer in a debugger,
/// like gdb.

/// So, append to types std::shared/weak/unique pointers and
/// functions std::make_shared/weak/unique the suffix _wrapper.
/// Do not use the prefix std::, of course. 

namespace  detail {


template<typename T>
struct pointer_target_wrapper
{
    using target_type = T;

    pointer_target_wrapper() : target{} {}
    template<typename... Args> pointer_target_wrapper(Args&&... args) : target(std::forward<Args>(args)...) {}

    target_type const* target_ptr() const { return &target; }
    target_type* target_ptr() { return &target; }

    target_type const* operator->() const { return &target; }
    target_type* operator->() { return &target; }
private:
    target_type target;
};


template<typename T, template<typename...> typename P>
struct smart_ptr_wrapper
{
    using target_type = T;
    using wrapper = pointer_target_wrapper<T>;
    using pointer_type = P<wrapper>;

    smart_ptr_wrapper() : pointer{} {}
    smart_ptr_wrapper(pointer_type  ptr) : pointer{} { pointer.swap(ptr); }
    template<typename Q> smart_ptr_wrapper(Q  ptr) : pointer{ ptr.wrapped_pointer() } {}

    smart_ptr_wrapper(std::nullptr_t) : pointer{nullptr} {}
    smart_ptr_wrapper&  operator=(std::nullptr_t) { pointer = nullptr; return *this; }
    bool  operator==(std::nullptr_t) const { return pointer == nullptr; }
    bool  operator!=(std::nullptr_t) const { return pointer != nullptr; }

    target_type const*  get() const { return pointer->target_ptr(); }
    target_type*  get() { return pointer->target_ptr(); }

    target_type const*  operator->() const { return get(); }
    target_type*  operator->() { return get(); }

    target_type const&  operator*() const { return *get(); }
    target_type&  operator*() { return *get(); }

    pointer_type const&  wrapped_pointer() const { return pointer; }
    pointer_type&  wrapped_pointer() { return pointer; }
private:
    pointer_type  pointer;
};


template<typename T> using shared_ptr_wrapper = smart_ptr_wrapper<T, std::shared_ptr>;
template<typename T> using weak_ptr_wrapper = smart_ptr_wrapper<T, std::weak_ptr>;
template<typename T> using unique_ptr_wrapper = smart_ptr_wrapper<T, std::unique_ptr>;


template<typename T, typename... Args>
shared_ptr_wrapper<T>  make_shared_wrapper(Args&&... args)
{
    return shared_ptr_wrapper<T>{ std::make_shared<typename shared_ptr_wrapper<T>::wrapper>(std::forward<Args>(args)...) };
}

template<typename T, typename... Args>
unique_ptr_wrapper<T>  make_unique_wrapper(Args&&... args)
{
    return unique_ptr_wrapper<T>{ std::make_unique<typename unique_ptr_wrapper<T>::wrapper>(std::forward<Args>(args)...) };
}


}

#   if PLATFORM() != PLATFORM_WINDOWS() && BUILD_RELEASE() == 0
#       define shared_ptr_wrapper detail::shared_ptr_wrapper
#       define weak_ptr_wrapper detail::weak_ptr_wrapper
#       define unique_ptr_wrapper detail::unique_ptr_wrapper
#       define make_shared_wrapper detail::make_shared_wrapper
#       define make_unique_wrapper detail::make_unique_wrapper
#   else
#       define shared_ptr_wrapper std::shared_ptr
#       define weak_ptr_wrapper std::weak_ptr
#       define unique_ptr_wrapper std::unique_ptr
#       define make_shared_wrapper std::make_shared
#       define make_unique_wrapper std::make_unique
#   endif

#endif
