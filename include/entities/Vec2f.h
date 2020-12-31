#ifndef FVM_VEC2F_H
#define FVM_VEC2F_H

#include <iostream>

namespace FVM { namespace Entities {

    /**
     * @brief The Vec2f class represents a 2-element vector.
     * 
     * This can be used for 2D coordinates, 2D directions, 2-component colours, etc.
     * Several arithmetical operations are defined for those vectors.
     */
    class Vec2f {
        public:
            double v[2]; /**< @brief Array of the 2 values in the vector.*/
        public:
            /**
             * @brief Construct a new Vec2f object with (0, 0).
             */
            Vec2f();

            /**
             * @brief Construct a new Vec2f object from 2 components.
             * 
             * @param x First component of the vector.
             * @param y Second component of the vector.
             */
            Vec2f(double x, double y); 

            /**
             * @brief Construct a new Vec2f object from one value.
             * 
             * @param x Value given to the two components of the vector.
             */
            explicit Vec2f(double x); 

            /**
             * @brief Accesses the selected component of the vector, returning a reference.
             * 
             * @param index Index of the component to access.
             * @return double& Reference to the selected component.
             */
            double &operator[](int index);

            /**
             * @brief Returns the selected component of the vector.
             * 
             * @param index Index of the component to return.
             * @return double Selected component.
             */
            double operator[](int index) const; 

            /**
             * @brief Accesses the selected component of the vector, returning a reference.
             * 
             * @param index Index of the component to access.
             * @return double& Reference to the selected component.
             */
            double &operator[](unsigned int index);

            /**
             * @brief Returns the selected component of the vector.
             * 
             * @param index Index of the component to return.
             * @return double Selected component.
             */
            double operator[](unsigned int index) const; 

            /**
             * @brief Multiplies all components of the vector by a factor.
             * 
             * Returns (x1*a, y1*a).
             * 
             * @param scale Factor used to multiply all components of the vector.
             * @return Vec2f Resulting vector, (x1*a, y1*a).
             */
            Vec2f operator*(double scale) const;

            /**
             * @brief Element-wise multiplication of two vectors.
             * 
             * Returns (x1*x2, y1*y2).
             * 
             * @param other Vector used to multiply.
             * @return Vec2f Resulting vector, (x1*x2, y1*y2).
             */
            Vec2f operator*(const Vec2f &other) const;

            /**
             * @brief Divides all components of the vector by a factor.
             * 
             * Returns (x1/a, y1/a).
             * 
             * @param scale Factor used to divide all components of the vector.
             * @return Vec2f Resulting vector, (x1/a, y1/a).
             */
            Vec2f operator/(double scale) const;

            /**
             * @brief Elements-wise division by the provided vector.
             * 
             * Returns (x1/x2, y1/y2).
             * 
             * @param other Vector used to divide the components of this vector.
             * @return Vec2f Resulting vector, (x1/x2, y1/y2).
             */
            Vec2f operator/(const Vec2f &other) const;

            /**
             * @brief Divides all components of the vector by a factor.
             * 
             * Returns (x1/a, y1/a).
             * 
             * @param scale Factor used to divide all components of the vector.
             * @return Vec2f Resulting vector, (x1/a, y1/a).
             */
            Vec2f operator/(size_t scale) const;

            /**
             * @brief Adds two vectors.
             * 
             * Returns (x1+x2, y1+y2).
             * 
             * @param other Vector added to this vector.
             * @return Vec2f Resulting vector, (x1+x2, y1+y2).
             */
            Vec2f operator+(const Vec2f &other) const;

            /**
             * @brief Adds a factor to all components of the vector.
             * 
             * Returns (x1+a, y1+a).
             * 
             * @param factor Factor added to all components of the vector.
             * @return Vec2f Resulting vector, (x1+a, y1+a).
             */
            Vec2f operator+(double factor) const;

            /**
             * @brief Substracts a vector from this vector.
             * 
             * Returns (x1-x2, y1-y2).
             * 
             * @param other Vector to substract from this vector.
             * @return Vec2f Resulting vector, (x1-x2, y1-y2).
             */
            Vec2f operator-(const Vec2f &other) const;

            /**
             * @brief Substracts a factor from all components of the vector.
             * 
             * Returns (x1-a, y1-a).
             * 
             * @param factor Factor substracted from all components of the vector.
             * @return Vec2f Resulting vector, (x1-a, y1-a).
             */
            Vec2f operator-(double factor) const;

            /**
             * @brief Returns the vector negated.
             * 
             * Returns (-x1, -y1).
             * 
             * @return Vec2f Resulting vector, (-x1, -y1).
             */
            Vec2f operator-() const; 

            /**
             * @brief Tests equality between two vectors.
             * 
             * @param other Vector used to test equality.
             * @return true All two components of the vectors are equal.
             * @return false At least one component of the vectors is different.
             */
            bool operator==(const Vec2f &other) const;

            /**
             * @brief In-place multiplies all components of the vector by a factor.
             * 
             * Becomes (x1*a, y1*a).
             * 
             * @param scale Factor used to multiply all components of the vector.
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &operator*=(double scale);

            /**
             * @brief In-place element-wise multiplication of the vector by another vector.
             * 
             * Becomes (x1*x2, y1*y2).
             * 
             * @param other Vector used to multiply.
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &operator*=(const Vec2f &other);

            /**
             * @brief In-place divides all components of the vector by a factor.
             * 
             * Becomes (x1/a, y1/a).
             * 
             * @param scale Factor used to divide all components of the vector.
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &operator/=(double scale);

            /**
             * @brief In-place divides all components of the vector by a factor.
             * 
             * Becomes (x1/a, y1/a).
             * 
             * @param scale Factor used to divide all components of the vector.
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &operator/=(size_t scale);

            /**
             * @brief In-place elements-wise division by the provided vector.
             * 
             * Becomes (x1/x2, y1/y2).
             * 
             * @param other Vector used to divide the components of this vector.
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &operator/=(const Vec2f &other);

            /**
             * @brief In-place addition of another vector.
             * 
             * Becomes (x1+x2, y1+y2).
             * 
             * @param other Vector added to this vector.
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &operator+=(const Vec2f &other);

            /**
             * @brief In-place adds a factor to all components of the vector.
             * 
             * Becomes (x1+a, y1+a).
             * 
             * @param factor Factor added to all components of the vector.
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &operator+=(double factor);

            /**
             * @brief In-place substracts a vector from this vector.
             * 
             * Becomes (x1-x2, y1-y2).
             * 
             * @param other Vector to substract from this vector.
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &operator-=(const Vec2f &other); 

            /**
             * @brief In-place substracts a factor from all components of the vector.
             * 
             * Becomes (x1-a, y1-a).
             * 
             * @param factor Factor substracted from all components of the vector.
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &operator-=(double factor);

            /**
             * @brief Sets the components of the vector to the minimum of its components and the other's.
             * 
             * Becomes (min(x1, x2), min(y1, y2))
             * 
             * @param other Vector to calculate minimum components with.
             * @return Vec2f& Reference to the vector, used to chain operations.
             */
            Vec2f &min(const Vec2f &other);

            /**
             * @brief Sets the components of the vector to the minimum of its components and the provided factor.
             * 
             * Becomes (min(x1, a), min(y1, a))
             * 
             * @param other Factor to calculate minimum with.
             * @return Vec2f& Reference to the vector, used to chain operations.
             */
            Vec2f &min(double other);

            /**
             * @brief Sets the components of the vector to the maximum of its components and the other's.
             * 
             * Becomes (max(x1, x2), max(y1, y2))
             * 
             * @param other Vector to calculate maximum components with.
             * @return Vec2f& Reference to the vector, used to chain operations.
             */
            Vec2f &max(const Vec2f &other);

            /**
             * @brief Sets the components of the vector to the maximum of its components and the provided factor.
             * 
             * Becomes (max(x1, a), max(y1, a))
             * 
             * @param other Factor to calculate maximum with.
             * @return Vec2f& Reference to the vector, used to chain operations.
             */
            Vec2f &max(double other);

            /**
             * @brief Returns a vector with the minimum components of this vector and another.
             * 
             * Returns (min(x1, x2), min(y1, y2))
             * 
             * @param other Vector to calculate minimum components with.
             * @return Vec2f Reference to the vector, used to chain operations.
             */
            Vec2f getMin(const Vec2f &other) const;

            /**
             * @brief Returns a vector with the minimum components of this vector and a factor.
             * 
             * Returns (min(x1, a), min(y1, a))
             * 
             * @param other Factor to calculate minimum with.
             * @return Vec2f Reference to the vector, used to chain operations.
             */
            Vec2f getMin(double other) const;

            /**
             * @brief Returns a vector with the maximum components of this vector and another.
             * 
             * Returns (max(x1, x2), max(y1, y2))
             * 
             * @param other Vector to calculate maximum components with.
             * @return Vec2f Reference to the vector, used to chain operations.
             */
            Vec2f getMax(const Vec2f &other) const;

            /**
             * @brief Returns a vector with the maximum components of this vector and a factor.
             * 
             * Returns (max(x1, a), max(y1, a))
             * 
             * @param other Factor to calculate maximum with.
             * @return Vec2f Reference to the vector, used to chain operations.
             */
            Vec2f getMax(double other) const;

            /**
             * @brief Returns the magnitude of the vector.
             * 
             * Returns the L2 norm of the vector: sqrt(x^2 + y^2).
             * 
             * @return double Magnitude of the vector.
             */
            double magnitude() const;

            /**
             * @brief Returns the squared magnitude of the vector.
             * 
             * Returns x^2 + y^2. Useful because it is much faster than the norm,
             * and can be used instead of it in some situations.
             * 
             * @return double Squared magnitude of the norm.
             */
            double magnitudeSquared() const;

            /**
             * @brief Returns a the normalized vector.
             * 
             * Divides all components of the vector by its magnitude.
             * 
             * @return Vec2f Normalized vector.
             */
            Vec2f normalize() const;

            /**
             * @brief Normalizes the vector in-place, dividing it by its norm.
             * 
             * Divides all components of the vector by its magnitude.
             * 
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &normalize_inplace();

            /**
             * @brief Computes the dot product of this vector and another.
             * 
             * Returns v1.v2
             * 
             * @param other Vector to dot with this one.
             * @return double Dot product of the two vectors.
             */
            double dot(const Vec2f &other) const;

            /**
             * @brief Changes the vector in-place to polar coordinates.
             * 
             * Assumes the vector is in cartesian coordinates.
             * 
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &to_sph();

            /**
             * @brief Changes the vector in-place to cartesian coordinates.
             * 
             * Assumes the vector is in polar coordinates.
             * 
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &to_xy();

            /**
             * @brief Changes the vector in-place to cartesian coordinates with arbitrary axises.
             * 
             * Assumes the vector is in polar coordinates.
             * 
             * @param ref1 Axis used for x.
             * @param ref2 Axis used for y.
             * @return const Vec2f& Reference to the vector, used to chain operations.
             */
            const Vec2f &to_xy_offset(const Vec2f &ref1, const Vec2f &ref2);

            /**
             * @brief Returns the vector in polar coordinates.
             * 
             * Assumes the vector is in cartesian coordinates.
             * 
             * @return Vec2f polar coordinates of the vector.
             */
            Vec2f get_sph() const;

            /**
             * @brief Returns the vector in cartesian coordinates.
             * 
             * Assumes the vector is in polar coordinates.
             * 
             * @return Vec2f Cartesian coordinates of the vector.
             */
            Vec2f get_xy() const;

            /**
             * @brief Returns the vector in cartesian coordinates with arbitrary axises.
             * 
             * Assumes the vector is in polar coordinates.
             * 
             * @param ref1 Axis used for x.
             * @param ref2 Axis used for y.
             * @return Vec2f Cartesian coordinates of the vector.
             */
            Vec2f get_xy_offset(const Vec2f &ref1, const Vec2f &ref2) const;

            /**
             * @brief Returns a vector of the natural logarithm of all components of the vector.
             * 
             * Returns (ln(x), ln(y))
             * 
             * @return Vec2f Vector made of the natural logarithm of all components of the vector.
             */
            Vec2f ln() const;

            /**
             * @brief Returns a vector of the square root of all components of the vector.
             * 
             * Returns (sqrt(x), sqrt(y))
             * 
             * @return Vec2f Vector made of the square root of all components of the vector.
             */
            Vec2f sqrt() const;

            /**
             * @brief Returns a vector of the exponential of all components of the vector.
             * 
             * Returns (e^x, e^y).
             * 
             * @return Vec2f Vector made of the exponential of all components of the vector.
             */
            Vec2f exp() const;

            /**
             * @brief Returns a vector of the components of the vector to the specified power.
             * 
             * Returns (x^a, y^a).
             * 
             * @param exp Power to be applied to all components.
             * @return Vec2f Vector made of the components of the vector to the specified power.
             */
            Vec2f pow(double exp) const;

            /**
             * @brief In-place raise the components of the vector to the specified power.
             * 
             * Becomes (x^a, y^a).
             * 
             * @param exp Power to be applied to all components.
             * @return Vec2f& Reference to the vector, used to chain operations.
             */
            Vec2f &pow_inplace(double exp);

            /**
             * @brief Returns a vector of the components of the vector rounded down.
             * 
             * Returns (floor(x), floor(y))
             * 
             * @return Vec2f Vector made of the components of the vector rounded down.
             */
            Vec2f floor() const;

            /**
             * @brief Returns a vector of the components of the vector rounded up.
             * 
             * Returns (ceil(x), ceil(y))
             * 
             * @return Vec2f Vector made of the components of the vector rounded up.
             */
            Vec2f ceil() const;

            /**
             * @brief In-place rounds the components to the nearest integer value.
             * 
             * Becomes (round(x), round(y))
             * 
             * @return Vec2f& Reference to the vector, used to chain operations.
             */
            Vec2f &round_inplace();

            /**
             * @brief In-place limits the components of the vector to a minimum and maximum value.
             * 
             * @param minimum Minimum value of the components.
             * @param maximum Maximum value of the components.
             * @return Vec2f& Reference to the vector, used to chain operations.
             */
            Vec2f &clamp(double minimum, double maximum);

            /**
             * @brief Returns the x component of the vector
             * 
             * @return double x component of the vector.
             */
            double x() const;
            
            /**
             * @brief Returns the y component of the vector
             * 
             * @return double y component of the vector.
             */
            double y() const;

            /**
             * @brief Returns a reference to the x component of the vector
             * 
             * @return double Reference to the x component of the vector.
             */
            double& x();
            
            /**
             * @brief Returns a reference to the y component of the vector
             * 
             * @return double Reference to the y component of the vector.
             */
            double& y();
    };
}}

/**
 * @brief Formats a vector to be displayed.
 * 
 * @param output Output stream.
 * @param v Vector to be displayed.
 * @return std::ostream& Output stream.
 */
std::ostream &operator<<(std::ostream &output, const FVM::Entities::Vec2f &v);

/**
 * @brief Multiplies a factor with a vector.
 * 
 * Returns (a*x, a*y).
 * 
 * @param factor Factor multiplying the vector.
 * @param v Vector to be multiplied.
 * @return FVM::Entities::Vec2f Resulting Vector, (a*x, a*y).
 */
FVM::Entities::Vec2f operator*(const double factor, const FVM::Entities::Vec2f &v);

#endif