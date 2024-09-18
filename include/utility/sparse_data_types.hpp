#ifndef UTILITY_SPARSE_DATA_TYPES_HPP_INCLUDED
#   define UTILITY_SPARSE_DATA_TYPES_HPP_INCLUDED

#   include <utility/basic_numeric_types.hpp>
#   include <unordered_map>


using scalar = float_64_bit;


struct  sparse_vector
{
    static inline sparse_vector  axis(std::size_t const  idx){ sparse_vector v; v[idx] = 1.0; return v; }

    scalar  operator()(std::size_t  idx) const;
    scalar&  operator[](std::size_t  idx);

    void  erase(std::size_t const  idx) { coords_.erase(idx); }

    sparse_vector&  operator+=(sparse_vector const&  v);
    sparse_vector&  operator-=(sparse_vector const&  v);
    sparse_vector&  operator*=(scalar  c);
    sparse_vector&  operator/=(scalar  c) { return this->operator*=(1.0 / c); }

    std::unordered_map<std::size_t, scalar> const&  coords() const { return coords_; }

private:
    std::unordered_map<std::size_t, scalar>  coords_{};
};

sparse_vector  operator-(sparse_vector const&  v);

sparse_vector  operator+(sparse_vector const&  u, sparse_vector const&  v);
sparse_vector  operator-(sparse_vector const&  u, sparse_vector const&  v);

scalar  dot(sparse_vector const&  u, sparse_vector const&  v);
scalar  length(sparse_vector const&  v);
sparse_vector  normalized(sparse_vector const&  v);

sparse_vector  operator*(scalar c, sparse_vector const&  v);
inline sparse_vector  operator*(sparse_vector const&  v, scalar const  c) { return c * v; }
inline sparse_vector  operator/(sparse_vector const&  v, scalar const  c) { return (1.0/c) * v; }


struct  sparse_orthogonal_basis
{
    explicit sparse_orthogonal_basis(std::size_t  num_dimensions) : dim_{ num_dimensions }, vectors_{} {}

    std::size_t  size() const { return dim_; }

    sparse_vector  operator()(std::size_t  idx) const;
    sparse_vector&  operator[](std::size_t  idx);

    sparse_vector  out(sparse_vector const& v) const;
    sparse_vector  in(sparse_vector const& v) const;

    void  erase(std::size_t const  idx) { vectors_.erase(idx); }

    std::unordered_map<std::size_t, sparse_vector> const&  vectors() const { return vectors_; }

private:
    std::size_t  dim_;
    std::unordered_map<std::size_t, sparse_vector>  vectors_{};
};


#endif
