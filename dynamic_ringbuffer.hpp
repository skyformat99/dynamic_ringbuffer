//
// dsa is a utility library of data structures and algorithms built with C++11.
// This file (dynamic_ringbuffer.hpp) is part of the dsa project.
//
// dynamic ringbuffer; a resizable, allocator-aware implementation of an
// STL-style circular buffer for C++11 and later.
//
// A description of the circular buffer data structure can be found here:
//
//      https://en.wikipedia.org/wiki/Circular_buffer
//
// author: Dalton Woodard
// contact: daltonmwoodard@gmail.com
// repository: https://github.com/daltonwoodard/dynamic_ringbuffer.git
// license:
//
// Copyright (c) 2016 DaltonWoodard. See the COPYRIGHT.md file at the top-level
// directory or at the listed source repository for details.
//
//      Licensed under the Apache License. Version 2.0:
//          https://www.apache.org/licenses/LICENSE-2.0
//      or the MIT License:
//          https://opensource.org/licenses/MIT
//      at the licensee's option. This file may not be copied, modified, or
//      distributed except according to those terms.
//

#ifndef DSA_DYNAMIC_RINGBUFFER_HPP
#define DSA_DYNAMIC_RINGBUFFER_HPP

#include <exception>    // std::runtime_error
#include <limits>       // std::numeric_limits
#include <memory>       // std::unique_ptr, std::allocator,
                        // std::allocator_traits
#include <stdexcept>    // std::length_error
#include <type_traits>  // std::remove_cv, std::is_nothrow_move_constructible,
                        // std::is_nothrow_copy_constructible,
                        // std::is_destructible,
                        // std::is_nothrow_destructible
#include <utility>      // std::forward, std::move, std::swap


namespace dsa
{
namespace
{
    template <typename T>
    struct memblock
    {
        alignas (alignof (T)) unsigned char data [sizeof (T)];
    };
}   // annonymous namespace

    /*
     *  Description
     *  -----------
     *
     *  dsa::dynamic_ringbuffer <> is an implementation of a resizable,
     *  allocator-aware STL-style circular buffer.
     *
     *  A description of the circular buffer data structure can be found here:
     *
     *      https://en.wikipedia.org/wiki/Circular_buffer
     *
     *  The size of the buffer is dynamic. Please see
     *  https://github.com/daltonwoodard/ringbuffer.git for a fixed-size
     *  implementation.
     *
     *  This implementation is NOT threadsafe. Please see
     *  https://github.com/daltonwoodard/atomic_dynamic_ringbuffer.git for a
     *  thread-safe version implemented with C++ atomic primitives.
     *
     *  If and only if the stored object type T has strong exception safe
     *  constructors is the following guaranteed:
     *
     *  If this object is successfully constructed, then throughout the
     *  object's lifetime:
     *
     *      i.  all view operations are guaranteed as nothrow
     *      ii. write operations provide the strong exception safety guarantee
     *
     *  Template Parameters
     *  -------------------
     *  - T: the object type to be buffered. This type does *not* have to be
     *  default constructible.
     *
     *  - Alloc: the allocator type used in internal buffer allocations. The
     *  default is std::allocator <T>. Allocators of type Alloc must:
     *
     *      i.  allocate blocks of memory of size:      sizeof (T)
     *      ii. allocate blocks of memory of alignment: alignof (T)
     *
     *  Class Scoped Enumerations
     *  -------------------------
     *  - dynamic_ringbuffer::resize_policy [default: resize]
     *      controls behavior of the container when capacity == 0: if the
     *      value is equal to no_resize then the behavior depends on the
     *      overwrite_policy; if the value is resize then the internal buffer
     *      is reallocated with a larger maximum capacity.
     *
     *  - dynamic_ringbuffer::overwrite_policy [default: no_overwrite]
     *      controls behavior of the container when capacity == 0 and the
     *      resize policy is no_resize: if the value is equal to overwrite,
     *      then upon a call to push/push_back or emplace/emplace_back the
     *      value held previously at the front of the buffer is overwritten;
     *      if the value is equal to no_overwrite then upon a call to
     *      push/push_back or emplace/emplace_back an exception of type
     *      std::runtime_error is emitted.
     *
     *  Member Types
     *  ------------
     *  - allocator_type:  Alloc;
     *  - value_type:      T;
     *  - size_type:       std::size_t;
     *  - difference_type: std::ptrdiff_t;
     *  - pointer:         value_type *;
     *  - const_pointer:   value_type const *;
     *  - reference:       value_type &;
     *  - const_reference: value_type const &;
     *
     *  - iterator:               models RandomAccessIterator
     *  - const_iterator:         models RandomAccessIterator
     *  - reverse_iterator:       std::reverse_iterator <iterator>;
     *  - const_reverse_iterator: std::reverse_iterator <const_iterator>;
     *
     *  Constructors
     *  ------------
     *  dynamic_ringbuffer (void):
     *      - default constructs an empty buffer; default constructs allocator;
     *        defaults resize and overwrite policies
     *      - nothrow if Alloc is nothrow default constructible
     *
     *  dynamic_ringbuffer (allocator_type const &):
     *      - default constructs an empty buffer; copy constructs allocator;
     *        defaults resize and overwrite policies
     *
     *  dynamic_ringbuffer (size_type,
     *                      enum resize_policy = default,
     *                      enum overwrite_policy = default):
     *      - constructs a buffer with capacity sufficienct to hold at least
     *        the provided size buffered elements; sets resize and overwrite
     *        policies
     *
     *  dynamic_ringbuffer (size_type,
     *                      allocator_type const &,
     *                      enum resize_policy = default,
     *                      enum overwrite_policy = default):
     *      - constructs a buffer with capacity sufficienct to hold at least
     *        the provided size buffered elements; copy constructs allocator;
     *        sets resize and overwrite policies
     *
     *  dynamic_ringbuffer (dynamic_ringbuffer const &):
     *      - copy constructs buffer; copy constructs allocator; copies resize
     *        and overwrite policies
     *
     *  dynamic_ringbuffer (dynamic_ringbuffer const &, allocator_type const &):
     *      - copy constructs buffer; copy constructs allocator; copies resize
     *        and overwrite policies
     *
     *  dynamic_ringbuffer (dynamic_ringbuffer &&):
     *      - move constructs buffer; move constructs allocator; copies resize
     *        and overwrite policies
     *      - nothrow if Alloc is nothrow move constructible
     *
     *  dynamic_ringbuffer (dynamic_ringbuffer &&, allocator_type const &):
     *      - move constructs buffer; copy constructs allocator; copies resize
     *        and overwrite policies
     *
     *  Assignment Operators
     *  --------------------
     *  operator= (dynamic_ringbuffer const &):
     *      - copy assigns buffer; copy assigns allocator if Alloc is propagate
     *      on container copy assignment; copies resize and overwrite policies
     *
     *  operator= (dynamic_ringbuffer &&):
     *      - move assigns buffer; move assigns allocator if Alloc is propagate
     *      on container move assignment; copies resize and overwrite policies
     *      - nothrow if T is nothrow destructible and either:
     *          i.  Alloc is propagate on container move assignment and Alloc is
     *              nothrow move assignable;
     *          ii. or Alloc is not propagate on container move assignment and
     *              Alloc is always equal
     *
     *  Member Functions
     *  ----------------
     *  - front:      access the first element
     *  - back:       access the last element
     *
     *  - empty:    checks whether the buffer is empty
     *  - size:     returns the number of buffered elements
     *  - max_size: returns the maximum possible number of elements
     *  - reserve:  reserves storage large enough for the given capacity
     *  - capacity: returns the number of elements that can be held in
     *              currently allocated storage
     *  - resize:   changes the number of elements stored
     *  - shrink_to_fit: reduces memory usage by freeing unused memory
     *
     *  - get_allocator: returns the allocator associated with the container
     *
     *  - set_resize_policy: sets the resize policy for the container
     *  - get_resize_policy: returns the resize policy for the container
     *
     *  - set_overwrite_policy: sets the overwrite policy for the container
     *  - get_overwrite_policy: returns the overwrite policy for the container
     *
     *  - clear:                clears the contents of the buffer
     *  - push/push_back:       inserts an element at the end
     *  - emplace/emplace_back: constructs an element in-place at the end
     *  - pop/pop_front:        removes the first element
     *
     *  - swap: swaps the contents. Template typename T must be Swappable.
     */
    template <typename T, typename Alloc = std::allocator <memblock <T>>>
    class dynamic_ringbuffer
    {
        static_assert (
            sizeof (T) ==
                sizeof (typename std::allocator_traits <Alloc>::value_type),
            "allocator of type template parameter typename Alloc must allocate"
            " blocks of memory of size equal to size of template typename"
            " parameter T"
        );

        static_assert (
            alignof (T) ==
                alignof (typename std::allocator_traits <Alloc>::value_type),
            "allocator of type template parameter typename Alloc must allocate"
            " blocks of memory of alignment equal to alignment of template"
            " typename parameter T"
        );

    public:
        enum class resize_policy
        {
            resize = 0,
            no_resize
        };

        enum class overwrite_policy
        {
            no_overwrite = 0,
            overwrite
        };

        using allocator_type = Alloc;

    private:
        /* allocator-relevant types derived from template parameter Alloc */
        using alloc_traits          = std::allocator_traits <allocator_type>;
        using alloc_value_type      = typename alloc_traits::value_type;
        using alloc_pointer         = typename alloc_traits::pointer;
        using alloc_difference_type = typename alloc_traits::difference_type;
        using alloc_size_type       = typename alloc_traits::size_type;

        using alloc_propagate_on_container_copy_assignment =
            typename alloc_traits::propagate_on_container_copy_assignment;

        using alloc_propagate_on_container_move_assignment =
            typename alloc_traits::propagate_on_container_move_assignment;

        using alloc_propagate_on_container_swap =
            typename alloc_traits::propagate_on_container_swap;

        using alloc_is_always_equal = typename alloc_traits::is_always_equal;

        /* general types derived from template parameter Alloc  */
        using alloc_is_nothrow_default_constructible =
            typename std::is_nothrow_default_constructible <
                allocator_type
            >::type;

        using alloc_is_nothrow_copy_constructible =
            typename std::is_nothrow_copy_constructible <
                allocator_type
            >::type;

        using alloc_is_nothrow_move_constructible =
            typename std::is_nothrow_move_constructible <
                allocator_type
            >::type;

        using alloc_is_nothrow_copy_assignable =
            typename std::is_nothrow_copy_assignable <
                allocator_type
            >::type;

        using alloc_is_nothrow_move_assignable =
            typename std::is_nothrow_move_assignable <
                allocator_type
            >::type;

        static constexpr auto buffer_growth_rate = 1.5;

        static std::size_t bump_up (std::size_t n) noexcept
        {
            static constexpr std::size_t max_size {
                std::numeric_limits <std::size_t>::max () / sizeof (T)
            };

            static constexpr std::size_t threshold {
                static_cast <std::size_t> (
                    max_size / buffer_growth_rate
                )
            };

            std::size_t r = 2;

            while (r < n) {
                if (r > threshold) {
                    r = max_size;
                    break;
                } else {
                    r *= buffer_growth_rate;
                }
            }

            return r;
        }

        using backing_pointer         = memblock <T> *;
        using backing_const_pointer   = memblock <T> const *;
        using backing_reference       = memblock <T> &;
        using backing_const_reference = memblock <T> const &;

        allocator_type _alloc;

        /* policies determining behavior when capacity == 0 */
        enum resize_policy _rspolicy;
        enum overwrite_policy _owpolicy;

        /* number of buffered elements */
        std::size_t _buffered;

        /* maximum currently allocated capacity */
        std::size_t _capacity;

        /* our actual buffer */
        memblock <T> * _buffer;

        /* _first and _last point to the first and last buffer positions */
        backing_pointer _first;
        backing_pointer _last;

        template <typename U>
        class iterator_impl;

        /*
         * The iterators _tail and _head are privileged in
         * that the logical space they work in is the whole buffer. They are
         * later used to represent the bounding logical regions when creating
         * iterators to the buffer.
         */
        iterator_impl <T> _head;
        iterator_impl <T> _tail;

        template <typename U>
        class iterator_impl : public std::iterator <
            std::random_access_iterator_tag, U, std::ptrdiff_t, U *, U &
        >
        {
        private:
            using iter_type = std::iterator <
                std::random_access_iterator_tag, U, std::ptrdiff_t, U *, U &
            >;

            U * _iter;

            /*
             * these pointers bound the logical region with the current state of
             * the circular buffer.
             */
            U * _lfirst;
            U * _llast;

            /* these pointers bound the address space of the backing buffer */
            U * _rfirst;
            U * _rlast;

            /*
             * the real size of the buffer; i.e. 1 plus the distance between
             * _rfirst and _rlast.
             */
            std::size_t _bufsize;

            /* checks whether the logical region is actually contiguous */
            bool logical_region_is_contiguous (void) const noexcept
            {
                return _lfirst <= _llast;
            }

        public:
            using difference_type   = typename iter_type::difference_type;
            using value_type        = typename iter_type::value_type;
            using pointer           = typename iter_type::pointer;
            using reference         = typename iter_type::reference;
            using iterator_category = typename iter_type::iterator_category;

            iterator_impl (void) noexcept
                : _iter    {nullptr}
                , _lfirst  {nullptr}
                , _llast   {nullptr}
                , _rfirst  {nullptr}
                , _rlast   {nullptr}
                , _bufsize {0}
            {}

            iterator_impl (pointer p,
                           pointer logical_first,
                           pointer logical_last,
                           pointer real_first,
                           pointer real_last)
                noexcept
                : _iter    {p}
                , _lfirst  {logical_first}
                , _llast   {logical_last}
                , _rfirst  {real_first}
                , _rlast   {real_last}
                , _bufsize {static_cast <std::size_t> (real_last - real_first)}
            {}

            iterator_impl & operator= (iterator_impl const & other) noexcept
                = default;

            void swap (iterator_impl & other) noexcept
            {
                std::swap (this->_iter, other._iter);
                std::swap (this->_first, other._first);
                std::swap (this->_last, other._last);
            }

            iterator_impl & operator++ (void)
            {
                if (_iter == _llast) {
                    _iter = _lfirst;
                } else {
                    _iter += 1;
                }

                return *this;
            }

            iterator_impl & operator-- (void)
            {
                if (_iter == _lfirst) {
                    _iter = _llast;
                } else {
                    _iter -= 1;
                }

                return *this;
            }

            iterator_impl operator++ (int)
            {
                auto const tmp {*this};

                if (_iter == _llast) {
                    _iter = _lfirst;
                } else {
                    _iter += 1;
                }

                return tmp;
            }

            iterator_impl operator-- (int)
            {
                auto const tmp {*this};

                if (_iter == _lfirst) {
                    _iter = _llast;
                } else {
                    _iter -= 1;
                }

                return tmp;
            }

            /*
             * the parameter n is assumed to be valid; that is, such that
             * the resulting pointer after [_iter += n] remains logically
             * bound between _lfirst and _llast.
             */
            iterator_impl & operator+= (difference_type n)
            {
                /*
                 * two cases:
                 *  i.  _lfirst <= _llast in the address space (non-wraparound)
                 *  ii. _llast < _lfirst in the address space (wraparound)
                 *      a. this is above _lfirst
                 *      b. this is below _llast
                 */

                /* case i. */
                if (this->logical_region_is_contiguous ()) {
                    _iter += n;
                /* cast ii.a. */
                } else if (_lfirst <= _iter) {
                    /* stays in-bounds */
                    if (_iter + n <= _rlast) {
                        _iter += n;
                    /* overflow past-the-end (note: n > 0) */
                    } else {
                        _iter = _rfirst + (n - 1 - (_rlast - _iter));
                    }
                /* cast ii.b. */
                } else {
                    /* stays in-bounds */
                    if (_iter + n >= _rfirst) {
                        _iter += n;
                    /* underflows before-the-beginning (note: n < 0) */
                    } else {
                        auto const m {-n};
                        _iter = _rlast - (m - 1 - (_iter - _rfirst));
                    }
                }

                return *this;
            }

            iterator_impl & operator-= (difference_type n)
            {
                return this->operator+= (-n);
            }

            iterator_impl operator+ (difference_type n) const
            {
                auto tmp = *this;
                return (tmp += n); 
            }

            iterator_impl operator- (difference_type n) const
            {
                auto tmp = *this;
                return (tmp -= n);
            }

            difference_type operator- (iterator_impl const & rhs) const
            {
                /* normal configuration -- non-wraparound case */
                if (_first < _last) {
                    return this->_iter - rhs._iter;
                /*
                 * _last is behind _first (in the address space) --
                 * wraparound case, so the space from _last to _first
                 * is uninitialized.
                 */
                } else {
                    /*
                     * three cases:
                     *  i.    both iters are above _lfirst
                     *  ii.   both iters are below _llast
                     *  iii.  iters are split above _lfirst and below _llast
                     *      a. this is above _lfirst and rhs is below _llast
                     *      b. rhs is above _lfirst and this is below _llast
                     */
                    if (_lfirst <= this->_iter && _lfirst <= rhs._iter) {
                        return this->_iter - rhs._iter;
                    } else if (this->_iter <= _llast && rhs._iter <= _llast) {
                        return this->_iter - rhs._iter;
                    } else if (_lfirst <= this->_iter && rhs._iter <= _llast) {
                        return (this->_iter - rhs._iter) -
                            static_cast <difference_type> (_bufsize);
                    /* if (_lfirst <= rhs._iter && this->_iter <= _llast) */
                    } else {
                        return static_cast <difference_type> (_bufsize) -
                            (rhs._iter - this->_iter);
                    }
                }
            }

            bool operator== (iterator_impl const & rhs) const
            {
                return this->_iter == rhs._iter;
            }

            bool operator!= (iterator_impl const & rhs) const
            {
                return this->_iter != rhs._iter;
            }

            bool operator< (iterator_impl const & rhs) const
            {
                /* logical: return this->_iter < rhs._iter; */
                return rhs - *this > 0;
            }

            bool operator> (iterator_impl const & rhs) const
            {
                /* logical: return _iter > rhs._iter; */
                return *this - rhs > 0;
            }

            bool operator<= (iterator_impl const & rhs) const
            {
                /* logical: return _iter <= rhs._iter; */
                if (*this == rhs) {
                    return true;
                } else {
                    return *this < rhs;
                }
            }

            bool operator>= (iterator_impl const & rhs) const
            {
                /* logical: return _iter >= rhs._iter; */
                if (*this == rhs) {
                    return true;
                } else {
                    return *this > rhs;
                }
            }

            reference operator* (void) const
            {
                return *_iter;
            }

            reference operator[] (difference_type n) const
            {
                return *(*this + n);
            }

            pointer addressof (void) const
            {
                return _iter;
            }
        };

        void set_buffer_pointers (void) noexcept
        {
            this->_first = &this->_buffer [0];
            this->_last = &this->_buffer [
                this->_capacity ? this->_capacity - 1 : 0
            ];
        }

        void set_buffer_iterators (void) noexcept
        {
            this->_head = iterator {
                reinterpret_cast <pointer> (this->_buffer),
                reinterpret_cast <pointer> (this->_first),
                reinterpret_cast <pointer> (this->_last),
                reinterpret_cast <pointer> (this->_first),
                reinterpret_cast <pointer> (this->_last)
            };
            this->_tail = iterator {
                reinterpret_cast <pointer> (this->_buffer),
                reinterpret_cast <pointer> (this->_first),
                reinterpret_cast <pointer> (this->_last),
                reinterpret_cast <pointer> (this->_first),
                reinterpret_cast <pointer> (this->_last)
            };
        }

        void reset (void) noexcept
        {
            this->_buffered = 0;
            this->_capacity = 0;
            this->_buffer = nullptr;
            this->_first = nullptr;
            this->_last = nullptr;
            this->_head = iterator {
                nullptr, nullptr, nullptr, nullptr, nullptr
            };
            this->_tail = iterator {
                nullptr, nullptr, nullptr, nullptr, nullptr
            };
        }

    public:
        using value_type      = T;
        using size_type       = std::size_t;
        using difference_type = std::ptrdiff_t;
        using pointer         = value_type *;
        using const_pointer   = value_type const *;
        using reference       = value_type &;
        using const_reference = value_type const &;

        using iterator        = iterator_impl <value_type>;
        using const_iterator  = iterator_impl <value_type const>;
        using reverse_iterator       = std::reverse_iterator <iterator>;
        using const_reverse_iterator = std::reverse_iterator <const_iterator>;

        dynamic_ringbuffer (void)
            noexcept (alloc_is_nothrow_default_constructible::value)
            : _alloc    {}
            , _rspolicy {resize_policy::resize}
            , _owpolicy {overwrite_policy::no_overwrite}
            , _buffered {0}
            , _capacity {0}
            , _buffer   {nullptr}
            , _first    {nullptr}
            , _last     {nullptr}
            , _head {}
            , _tail {}
        {}

        explicit dynamic_ringbuffer (allocator_type const & alloc)
            noexcept (alloc_is_nothrow_copy_constructible::value)
            : _alloc    {alloc}
            , _rspolicy {resize_policy::resize}
            , _owpolicy {overwrite_policy::no_overwrite}
            , _buffered {0}
            , _capacity {0}
            , _buffer   {nullptr}
            , _first    {nullptr}
            , _last     {nullptr}
            , _head  {}
            , _tail  {}
        {}

        dynamic_ringbuffer (
            size_type count,
            enum resize_policy rspol = resize_policy::resize,
            enum overwrite_policy owpol = overwrite_policy::no_overwrite
        )
            : _alloc    {}
            , _rspolicy {rspol}
            , _owpolicy {owpol}
            , _buffered {0}
            , _capacity {bump_up (count)}
            , _buffer   {alloc_traits::allocate (_alloc, _capacity)}
        {
            this->set_buffer_pointers ();
            this->set_buffer_iterators ();
        }

        dynamic_ringbuffer (
            size_type count,
            allocator_type const & alloc,
            enum resize_policy rspol = resize_policy {},
            enum overwrite_policy owpol = overwrite_policy {}
        )
            : _alloc    {alloc}
            , _rspolicy {rspol}
            , _owpolicy {owpol}
            , _buffered {0}
            , _capacity {bump_up (count)}
            , _buffer   {alloc_traits::allocate (_alloc, _capacity)}
        {
            this->set_buffer_pointers ();
            this->set_buffer_iterators ();
        }

        dynamic_ringbuffer (dynamic_ringbuffer const & other)
            : _alloc {
                alloc_traits::select_on_container_copy_construction (
                    other._alloc
                )
            }
            , _rspolicy {other._rspolicy}
            , _owpolicy {other._owpolicy}
            , _buffered {other._buffered}
            , _capacity {other._capacity}
            , _buffer   {alloc_traits::allocate (_alloc, _capacity)}
        {
            this->set_buffer_pointers ();
            this->set_buffer_iterators ();

            auto ti {this->_tail};
            auto oi {other.cbegin ()};
            auto ob {other._buffered};

            while (ob) {
                auto const addr {ti.addressof ()};
                new (addr) value_type {*oi};
                ti += 1;
                oi += 1;
                ob -= 1;
            }

            this->_tail += this->_buffered;
        }

        dynamic_ringbuffer (
            dynamic_ringbuffer const & other, allocator_type const & alloc
        )
            : _alloc    {alloc}
            , _rspolicy {other._rspolicy}
            , _owpolicy {other._owpolicy}
            , _buffered {other._buffered}
            , _capacity {other._capacity}
            , _buffer   {alloc_traits::allocate (_alloc, _capacity)}
        {
            this->set_buffer_pointers ();
            this->set_buffer_iterators ();

            auto ti {this->_tail};
            auto oi {other.cbegin ()};
            auto ob {other._buffered};

            while (ob) {
                auto const addr {ti.addressof ()};
                new (addr) value_type {*oi};
                ti += 1;
                oi += 1;
                ob -= 1;
            }

            this->_tail += this->_buffered;
        }

        dynamic_ringbuffer (dynamic_ringbuffer && other)
            noexcept (alloc_is_nothrow_move_constructible::value)
            : _alloc    {std::move (other._alloc)}
            , _rspolicy {other._rspolicy}
            , _owpolicy {other._owpolicy}
            , _buffered {other._buffered}
            , _capacity {other._capacity}
            , _buffer   {other._buffer}
        {
            this->set_buffer_pointers ();
            this->set_buffer_iterators ();

            other.reset ();
        }

        dynamic_ringbuffer (
            dynamic_ringbuffer && other, allocator_type const & alloc
        )
            : _alloc    {alloc}
            , _rspolicy {other._rspolicy}
            , _owpolicy {other._owpolicy}
            , _buffered {other._buffered}
            , _capacity {other._capacity}
        {
            if (this->_alloc == other._alloc) {
                this->_buffer = other._buffer;
                this->_head = other._head;
                this->_tail = other._tail;
                this->set_buffer_iterators ();
            } else {
                this->_buffer = alloc_traits::allocate (
                    this->_alloc, this->_capacity
                );

                this->set_buffer_pointers ();
                this->set_buffer_iterators ();

                auto ti {this->_tail};
                auto oi {other.cbegin ()};
                auto ob {other._buffered};

                while (ob) {
                    auto const addr {ti.addressof ()};
                    new (addr) value_type {std::move (*oi)};
                    ti += 1;
                    oi += 1;
                    ob -= 1;
                }

                this->_tail += this->_buffered;
                other.clear ();
            }

            other.reset ();
        }

        dynamic_ringbuffer & operator= (dynamic_ringbuffer const & other)
        {
            this->clear ();
            this->_rspolicy = other._rspolicy;
            this->_owpolicy = other._owpolicy;

            if (alloc_propagate_on_container_copy_assignment::value) {
                if (this->_alloc != other._alloc && this->_buffer != nullptr) {
                    alloc_traits::deallocate (
                        this->_alloc,
                        this->_buffer,
                        this->_capacity
                    );
                }

                this->_alloc = other._alloc;
                this->_capacity = other._capacity;
                this->_buffer = alloc_traits::allocate (
                    this->_alloc, this->_capacity
                );
            } else {
                if (this->_capacity != other._capacity) {
                    if (this->_buffer != nullptr) {
                        alloc_traits::deallocate (
                            this->_alloc,
                            this->_buffer,
                            this->_capacity
                        );
                    }

                    this->_capacity = other._capacity;
                    this->_buffer = alloc_traits::allocate (
                        this->_alloc, this->_capacity
                    );
                }
            }

            this->set_buffer_pointers ();
            this->set_buffer_iterators ();

            auto oi {other.cbegin ()};
            auto ob {other._buffered};

            while (ob) {
                auto const addr {this->_tail.addressof ()};
                new (addr) value_type {*oi};
                oi += 1;
                ob -= 1;

                /*
                 * we leave this here in the case that T is not nothrow copy
                 * constructible, so that in the case of an exception being
                 * thrown the container remains in a consistent state.
                 */
                this->_buffered += 1;
                this->_tail += 1;
            }

            return *this;
        }

        dynamic_ringbuffer & operator= (dynamic_ringbuffer && other)
            noexcept (
                std::is_nothrow_destructible <value_type>::value &&
                (
                    (alloc_propagate_on_container_move_assignment::value &&
                     alloc_is_nothrow_move_assignable::value) ||
                    (not alloc_propagate_on_container_move_assignment::value &&
                     alloc_is_always_equal::value)
                )
            )
        {
            this->clear ();
            this->_rspolicy = other._rspolicy;
            this->_owpolicy = other._owpolicy;

            if (alloc_propagate_on_container_move_assignment::value) {
                if (this->_buffer != nullptr) {
                    alloc_traits::deallocate (
                        this->_alloc,
                        this->_buffer,
                        this->_capacity
                    );
                }

                this->_alloc = std::move (other._alloc);
                this->_buffered = other._buffered;
                this->_capacity = other._capacity;
                this->_buffer = other._buffer;
                this->set_buffer_pointers ();
                this->set_buffer_iterators ();
            } else if (this->_alloc == other._alloc) {
                if (this->_buffer != nullptr) {
                    alloc_traits::deallocate (
                        this->_alloc,
                        this->_buffer,
                        this->_capacity
                    );
                }

                this->_buffered = other._buffered;
                this->_capacity = other._capacity;
                this->_buffer = other._buffer;
                this->set_buffer_pointers ();
                this->set_buffer_iterators ();
            } else {
                if (this->_capacity < other._buffered) {
                    if (this->_buffer != nullptr) {
                        alloc_traits::deallocate (
                            this->_alloc,
                            this->_buffer,
                            this->_capacity
                        );
                    }

                    this->_capacity = other._capacity;
                    this->_buffer = alloc_traits::allocate (
                        this->_alloc, this->_capacity
                    );
                    this->set_buffer_pointers ();
                    this->set_buffer_iterators ();
                }

                auto oi {other.cbegin ()};
                auto ob {other._buffered};

                while (ob) {
                    auto const addr {this->_tail.addressof ()};
                    new (addr) value_type {std::move (*oi)};
                    oi += 1;
                    ob -= 1;

                    /*
                     * we leave this here in the case that T is not nothrow move
                     * constructible, so that in the case of an exception being
                     * thrown the container remains in a consistent state.
                     */
                    this->_buffered += 1;
                    this->_tail += 1;
                }

                other.clear ();
                alloc_traits::deallocate (
                    other._alloc, other._buffer, other._capacity
                );
            }

            other.reset ();
            return *this;
        }

    private:
        static void destruct (pointer p)
            noexcept (std::is_nothrow_destructible <value_type>::value)
        {
            p->~value_type ();
        }

    public:
        ~dynamic_ringbuffer (void)
            noexcept (std::is_nothrow_destructible <value_type>::value)
        {
            if (this->_buffer) {
                auto it {this->end ()};

                while (_buffered > 0) {
                    it -= 1;
                    _buffered -= 1;
                    destruct (it.addressof ());
                }

                alloc_traits::deallocate (
                    this->_alloc, this->_buffer, this->_capacity
                );
            }
        }

        /* swaps the contents of the buffer */
        void swap (dynamic_ringbuffer & other)
            noexcept (
                not alloc_propagate_on_container_swap::value ||
                // TODO: replace this with a portable custom implementation
                std::__is_nothrow_swappable <allocator_type>::value
            )
        {
            if (alloc_propagate_on_container_swap::value) {
                using std::swap;
                swap (this->_alloc, other._alloc);
            }

            std::swap (this->_buffer, other._buffer);
            std::swap (this->_buffered, other._buffered);
            std::swap (this->_capacity, other._capacity);
            std::swap (this->_rspolicy, other._rspolicy);
            std::swap (this->_owpolicy, other._owpolicy);
            std::swap (this->_first, other._first);
            std::swap (this->_last, other._last);
            std::swap (this->_head, other._head);
            std::swap (this->_tail, other._tail);
        }

        /* returns the allocator associated with the container */
        allocator_type get_allocator (void) const
        {
            return _alloc;
        }

        /* checks whether the buffer is empty */
        bool empty (void) const noexcept
        {
            return _buffered != 0;
        }

        /* returns the number of elements stored in the buffer */
        std::size_t size (void) const noexcept
        {
            return _buffered;
        }

        /* returns the number of available spaces to add elements */
        std::size_t available (void) const noexcept
        {
            return _capacity - _buffered;
        }

        /*
         * returns the number of elements that can be held in the currently
         * allocated storage
         */
        std::size_t capacity (void) const noexcept
        {
            return _capacity;
        }

        /* returns the maximum number of elements that can be buffered */
        constexpr std::size_t max_size (void) const noexcept
        {
            return std::numeric_limits <std::size_t>::max () /
                sizeof (value_type);
        }

        /* reserves storage large enough to store at least new_cap elements */
        void reserve (std::size_t new_cap)
        {
            if (new_cap > _capacity) {
                if (new_cap > this->max_size ()) {
                    throw std::length_error {
                        "requested reserve of capacity larger than maximum size"
                    };
                } else {
                    auto const cap {bump_up (new_cap)};
                    auto const new_alloc {
                        alloc_traits::allocate (this->_alloc, cap)
                    };

                    {
                        auto insert_ptr {
                            reinterpret_cast <pointer> (new_alloc)
                        };
                        for (auto & e : *this) {
                            new (insert_ptr) value_type {std::move (e)};
                            insert_ptr += 1;
                        }
                    }

                    if (this->_buffer != nullptr) {
                        alloc_traits::deallocate (
                            this->_alloc,
                            this->_buffer,
                            this->_capacity
                        );
                    }

                    this->_capacity = cap;
                    this->_buffer = new_alloc;
                    this->set_buffer_pointers ();
                    this->set_buffer_iterators ();
                }
            }
        }

        /*
         * resizes the container to contain count elements
         *      if the current size is less than count, additional default
         *      inserted elements are appended.
         *
         *      if the current size is greater than count, the container is
         *      reduced to its first count elements.
         */
        void resize (std::size_t count)
        {
            if (count == 0) {
                this->clear ();
            } else {
                if (_capacity < count) {
                    this->reserve (count);
                }

                if (_buffered <= count) {
                    while (_buffered < count) {
                        auto const addr {_tail.addressof ()};
                        new (addr) value_type {};
                        _tail += 1;
                        _buffered += 1;
                    }
                } else {
                    while (_buffered > count) {
                        _tail -= 1;
                        _buffered -= 1;
                        destruct (_tail.addressof ());
                    }
                }
            }
        }

        /*
         * resizes the container to contain count elements:
         *      if the current size is less than count, additional elements
         *      are appended and initialized with copies of value;
         *
         *      if the current size is greater than count, the container is
         *      reduced to its first count elements.
         */
        void resize (std::size_t count, value_type const & value)
        {
            if (count == 0) {
                this->clear ();
            } else {
                if (_capacity < count) {
                    this->reserve (count);
                }

                if (_buffered <= count) {
                    while (_buffered < count) {
                        auto const addr {_tail.addressof ()};
                        new (addr) value_type {value};
                        _tail += 1;
                        _buffered += 1;
                    }
                } else {
                    while (_buffered > count) {
                        _tail -= 1;
                        _buffered -= 1;
                        destruct (_tail.addressof ());
                    }
                }
            }
        }

        /* requests the removal of unused capacity */
        void shrink_to_fit (void)
        {
            auto const bu {bump_up (_buffered)};
            if (bu < _capacity) {
                auto const new_alloc {
                    alloc_traits::allocate (this->_alloc, bu)
                };

                {
                    auto insert_ptr {reinterpret_cast <pointer> (new_alloc)};
                    for (auto & e : *this) {
                        new (insert_ptr) value_type {std::move (e)};
                        insert_ptr += 1;
                    }
                }

                if (this->_buffer != nullptr) {
                    alloc_traits::deallocate (
                        this->_alloc,
                        this->_buffer,
                        this->_capacity
                    );
                }

                this->_capacity = bu;
                this->_buffer = new_alloc;
                this->set_buffer_pointers ();
                this->set_buffer_iterators ();
            }
        }

        /* set the resize policy for the buffer */
        void set_resize_policy (enum resize_policy pol) noexcept
        {
            this->_rspolicy = pol;
        }

        /* returns the resize policy for the buffer */
        enum resize_policy get_resize_policy (void) const noexcept
        {
            return this->_rspolicy;
        }

        /* set the overwrite policy for the buffer */
        void set_overwrite_policy (enum overwrite_policy pol) noexcept
        {
            this->_owpolicy = pol;
        }

        /* returns the overwrite policy for the buffer */
        enum overwrite_policy get_overwrite_policy (void) const noexcept
        {
            return this->_owpolicy;
        }

        /* returns an iterator the start of the buffer */
        iterator begin (void) noexcept
        {
            return iterator {
                _head.addressof (),
                _head.addressof (),
                _tail == _head ? _tail.addressof ()
                               : (_tail - 1).addressof (),
                reinterpret_cast <pointer> (_first),
                reinterpret_cast <pointer> (_last)
            };
        }

        /* returns an iterator the end of the buffer */
        iterator end (void) noexcept
        {
            return iterator {
                _tail.addressof (),
                _head.addressof (),
                _tail == _head ? _tail.addressof ()
                               : (_tail - 1).addressof (),
                reinterpret_cast <pointer> (_first),
                reinterpret_cast <pointer> (_last)
            };
        }

        /* returns a const iterator the start of the buffer */
        const_iterator begin (void) const noexcept
        {
            return this->cbegin ();
        }

        /* returns a const iterator the end of the buffer */
        const_iterator end (void) const noexcept
        {
            return this->cend ();
        }

        /* returns a const iterator the start of the buffer */
        const_iterator cbegin (void) const noexcept
        {
            return const_iterator {
                _head.addressof (),
                _head.addressof (),
                _tail == _head ? _tail.addressof ()
                               : (_tail - 1).addressof (),
                reinterpret_cast <pointer> (_first),
                reinterpret_cast <pointer> (_last)
            };
        }

        /* returns a const iterator the end of the buffer */
        const_iterator cend (void) const noexcept
        {
            return const_iterator {
                _tail.addressof (),
                _head.addressof (),
                _tail == _head ? _tail.addressof ()
                               : (_tail - 1).addressof (),
                reinterpret_cast <pointer> (_first),
                reinterpret_cast <pointer> (_last)
            };
        }

        /* returns a reverse iterator to the start of the reversed buffer */
        reverse_iterator rbegin (void) noexcept
        {
            return reverse_iterator {this->end ()};
        }

        /* returns a reverse iterator to the end of the reversed buffer */
        reverse_iterator rend (void) noexcept
        {
            return reverse_iterator {this->begin ()};
        }

        /*
         * returns a const reverse iterator to the start of the reversed buffer
         */
        const_reverse_iterator rbegin (void) const noexcept
        {
            return this->crbegin ();
        }

        /* returns a const reverse iterator to the end of the reversed buffer */
        const_reverse_iterator rend (void) const noexcept
        {
            return this->crend ();
        }

        /*
         * returns a const reverse iterator to the start of the reversed buffer
         */
        const_iterator crbegin (void) const noexcept
        {
            return const_reverse_iterator {this->cend ()};
        }

        /* returns a const reverse iterator to the end of the reversed buffer */
        const_iterator crend (void) const noexcept
        {
            return const_reverse_iterator {this->cbegin ()};
        }

        /* returns a reference to the first element in the buffer */
        reference front (void) noexcept
        {
            return *_head;
        }

        /* returns a reference to the first element in the buffer */
        const_reference front (void) const noexcept
        {
            return *_head;
        }

        /* returns a reference to the last element in the buffer */
        reference back (void) noexcept
        {
            return _tail == _head ? _tail [0] : _tail [-1];
        }

        /* returns a reference to the last element in the buffer */
        const_reference back (void) const noexcept
        {
            return _tail == _head ? _tail [0] : _tail [-1];
        }

        /*
         * clears the contents of the buffer; the elements are guaranteed to be
         * destructed in the reverse-order they were added.
         */
        void clear (void)
            noexcept (std::is_nothrow_destructible <value_type>::value)
        {
            auto it {_tail};

            while (_buffered > 0) {
                it -= 1;
                _buffered -= 1;
                destruct (it.addressof ());
            }

            _tail = _head;
        }

        /*
         * Adds an object to the buffer if room is available.
         *
         * If no capacity is available, then:
         *      If the overwrite policy is set to no_overwrite this method
         *      throws an exception of type std::runtime_error.
         *
         *      If the overwrite policy is set to overwrite, this method
         *      overwrites the first element of the buffer.
         */
        void push (value_type const & v)
        {
            if (_buffered < _capacity) {
                auto const addr {_tail.addressof ()};
                new (addr) value_type {v};
                _tail += 1;
                _buffered += 1;
            } else {
                if (_rspolicy == resize_policy::resize) {
                    this->reserve (bump_up (_capacity));
                    auto const addr {_tail.addressof ()};
                    new (addr) value_type {v};
                    _tail += 1;
                    _buffered += 1;
                } else if (_owpolicy == overwrite_policy::overwrite) {
                    auto const addr {_tail.addressof ()};
                    destruct (addr);
                    new (addr) value_type {v};
                    _tail += 1;
                    _head += 1;
                } else {
                    throw std::runtime_error {"push back on full buffer"};
                }
            }
        }

        /*
         * Adds an object to the buffer if room is available.
         *
         * If no capacity is available, then:
         *      If the overwrite policy is set to no_overwrite this method
         *      throws an exception of type std::runtime_error.
         *
         *      If the overwrite policy is set to overwrite, this method
         *      overwrites the first element of the buffer.
         */
        void push (value_type && v)
        {
            if (_buffered < _capacity) {
                auto const addr {_tail.addressof ()};
                new (addr) value_type {std::move (v)};
                _tail += 1;
                _buffered += 1;
            } else {
                if (_rspolicy == resize_policy::resize) {
                    this->reserve (bump_up (_capacity));
                    auto const addr {_tail.addressof ()};
                    new (addr) value_type {std::move (v)};
                    _tail += 1;
                    _buffered += 1;
                } else if (_owpolicy == overwrite_policy::overwrite) {
                    auto const addr {_tail.addressof ()};
                    destruct (addr);
                    new (addr) value_type {std::move (v)};
                    _tail += 1;
                    _head += 1;
                } else {
                    throw std::runtime_error {"push back on full buffer"};
                }
            }
        }

        void push_back (value_type const & v)
        {
            this->push (v);
        }

        void push_back (value_type && v)
        {
            this->push (std::move (v));
        }

        /*
         * Adds an object to the buffer if room is available using in-place
         * construction.
         *
         * If no capacity is available, then:
         *      If the overwrite policy is set to no_overwrite this method
         *      throws an exception of type std::runtime_error.
         *
         *      If the overwrite policy is set to overwrite, this method
         *      overwrites the first element of the buffer.
         */
        template <typename ... Args>
        void emplace (Args && ... args)
        {
            if (_buffered < _capacity) {
                auto const addr {_tail.addressof ()};
                new (addr) value_type {std::forward <Args> (args)...};
                _tail += 1;
                _buffered += 1;
            } else {
                if (_rspolicy == resize_policy::resize) {
                    this->reserve (bump_up (_capacity));
                    auto const addr {_tail.addressof ()};
                    new (addr) value_type {std::forward <Args> (args)...};
                    _tail += 1;
                    _buffered += 1;
                } else if (_owpolicy == overwrite_policy::overwrite) {
                    auto const addr {_tail.addressof ()};
                    destruct (addr);
                    new (addr) value_type {std::forward <Args> (args)...};
                    _tail += 1;
                    _head += 1;
                } else {
                    throw std::runtime_error {"emplace on full buffer"};
                }
            }
        }

        template <typename ... Args>
        void emplace_back (Args && ... args)
        {
            this->emplace (std::forward <Args> (args)...);
        }

        /*
         * removes the first element from buffer if such an element exists, and
         * otherwise does nothing.
         */
        void pop (void)
            noexcept (std::is_nothrow_destructible <value_type>::value)
        {
            if (_buffered > 0) {
                destruct (_head.addressof ());
                _head += 1;
                _buffered -= 1;
            }
        }
    };
}   // namespace dsa

#endif // ifndef DSA_DYNAMIC_RINGBUFFER_HPP
