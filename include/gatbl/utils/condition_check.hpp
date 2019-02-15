#ifndef condition_check_HPP
#define condition_check_HPP

#include "gatbl/common.hpp"

namespace gatbl {
namespace utils {

template<bool debug = DEBUG>
struct condition_check;

template<>
struct condition_check<false>
{
    void set(bool) const {}
    void unchecked() const {}
    void checked() const {}
    void check() const {}
};

template<>
struct condition_check<true>
{
    void set(bool checked) const { _checked = checked; }
    void unchecked() const { _checked = false; }
    void checked() const { _checked = true; }
    void check(const char file[] = __FILE__, const char fun[] = __PRETTY_FUNCTION__, unsigned line = __LINE__) const
    {
        if (unlikely(!_checked)) {
            abort_message("%s:%s:%u: condition check failed", file, fun, line);
        }
    }

  private:
    mutable bool _checked = false;
};

} // namespace utils
} // namespace gatbl

#endif // condition_check_HPP
