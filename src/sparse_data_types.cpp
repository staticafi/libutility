#include <utility/sparse_data_types.hpp>
#include <utility/assumptions.hpp>
#include <utility/invariants.hpp>
#include <utility/development.hpp>
#include <unordered_set>
#include <algorithm>
#include <cmath>


scalar  sparse_vector::operator()(std::size_t const  idx) const
{
    auto const  it = coords_.find(idx);
    return it == coords_.end() ? 0.0 : it->second;
}


scalar&  sparse_vector::operator[](std::size_t const  idx)
{
    return coords_.insert({ idx, 0.0 }).first->second;
}


sparse_vector&  sparse_vector::operator+=(sparse_vector const&  v)
{
    for (auto  it = v.coords().begin(); it != v.coords().end(); ++it)
        this->operator[](it->first) += it->second;
    return *this;
}


sparse_vector&  sparse_vector::operator-=(sparse_vector const&  v)
{
    for (auto  it = v.coords().begin(); it != v.coords().end(); ++it)
        this->operator[](it->first) -= it->second;
    return *this;
}


sparse_vector&  sparse_vector::operator*=(scalar const  c)
{
    for (auto  it = coords_.begin(); it != coords_.end(); ++it)
        it->second *= c;
    return *this;
}


sparse_vector  operator-(sparse_vector const&  v)
{
    sparse_vector  w;
    for (auto  it = v.coords().begin(); it != v.coords().end(); ++it)
        w[it->first] = -it->second;
    return w;
}


sparse_vector  operator+(sparse_vector const&  u, sparse_vector const&  v)
{
    sparse_vector  w{ u };
    for (auto  it = v.coords().begin(); it != v.coords().end(); ++it)
        w[it->first] += it->second;
    return w;
}


sparse_vector  operator-(sparse_vector const&  u, sparse_vector const&  v)
{
    sparse_vector  w{ u };
    for (auto  it = v.coords().begin(); it != v.coords().end(); ++it)
        w[it->first] -= it->second;
    return w;
}


scalar  dot(sparse_vector const&  u, sparse_vector const&  v)
{
    scalar  w;
    for (auto  it = v.coords().begin(); it != v.coords().end(); ++it)
        w += u(it->first) * it->second;
    return w;
}


scalar  length(sparse_vector const&  v)
{
    return std::sqrt(dot(v,v));
}


sparse_vector  normalized(sparse_vector const&  v)
{
    return v / length(v);
}


sparse_vector  operator*(scalar c, sparse_vector const&  v)
{
    sparse_vector  w;
    for (auto  it = v.coords().begin(); it != v.coords().end(); ++it)
        w[it->first] = c * it->second;
    return w;
}


void  decompose(sparse_vector const&  v, sparse_vector const&  u, sparse_vector&  v_parallel_to_u, sparse_vector&  v_orthogonal_to_u)
{
    v_parallel_to_u = (dot(v, u) / dot(u, u)) * u;
    v_orthogonal_to_u = v - v_parallel_to_u;
}


bool  isfinite(sparse_vector const&  v)
{
    for (auto  it = v.coords().begin(); it != v.coords().end(); ++it)
        if (!std::isfinite(it->second))
            return false;
    return true;
}


sparse_vector  sparse_orthogonal_basis::operator()(std::size_t const  idx) const
{
    auto const  it = vectors_.find(idx);
    return it == vectors_.end() ? sparse_vector::axis(idx) : it->second;
}


sparse_vector&  sparse_orthogonal_basis::operator[](std::size_t const  idx)
{
    return vectors_.insert({ idx, sparse_vector{} }).first->second;
}


sparse_vector  sparse_orthogonal_basis::out(sparse_vector const& v) const
{
    sparse_vector  w;
    for (auto  it = v.coords().begin(); it != v.coords().end(); ++it)
        w += it->second * this->operator()(it->first);
    return w;
}


sparse_vector  sparse_orthogonal_basis::in(sparse_vector const& v) const
{
    sparse_vector  w;

    if (v.coords().empty())
        return w;

    std::unordered_set<std::size_t>  indices;
    for (auto  it = v.coords().begin(); it != v.coords().end(); ++it)
        indices.insert(it->first);
    for (auto  it = vectors_.begin(); it != vectors_.end(); ++it)
        indices.insert(it->first);

    for (std::size_t  i : indices)
    {
        sparse_vector const u{ this->operator()(i) };
        scalar const  c{ dot(u,v) / dot(u,u) };
        if (c != 0.0)
            w[i] = c;
    }
    return w;
}
