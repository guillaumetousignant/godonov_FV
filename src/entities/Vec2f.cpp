#include "entities/Vec2f.h"
#include <algorithm>
#include <cmath>

using FVM::Entities::Vec2f;

Vec2f::Vec2f() : v{0, 0} {
} 
Vec2f::Vec2f(double x, double y) : v{x, y} {
} 
Vec2f::Vec2f(double x) : v{x, x} {
} 
double &Vec2f::operator[](int index) {
    return v[index];
} 
double Vec2f::operator[](int index) const {
    return v[index];
} 
double &Vec2f::operator[](unsigned int index) {
    return v[index];
} 
double Vec2f::operator[](unsigned int index) const {
    return v[index];
} 
Vec2f Vec2f::operator*(double scale) const {
    return Vec2f(v[0] * scale, v[1] * scale);
} 
Vec2f Vec2f::operator*(const Vec2f &other) const {
    return Vec2f(v[0] * other[0], v[1] * other[1]);
} 
Vec2f Vec2f::operator/(double scale) const {
    return Vec2f(v[0] / scale, v[1] / scale);
} 
Vec2f Vec2f::operator/(size_t scale) const {
    return Vec2f(v[0] / scale, v[1] / scale);
} 
Vec2f Vec2f::operator/(const Vec2f &other) const {
    return Vec2f(v[0] / other[0], v[1] / other[1]);
}
Vec2f Vec2f::operator+(const Vec2f &other)  const{
    return Vec2f(v[0] + other.v[0], v[1] + other.v[1]);
} 
Vec2f Vec2f::operator+(double factor)  const{
    return Vec2f(v[0] + factor, v[1] + factor);
} 
Vec2f Vec2f::operator-(const Vec2f &other) const {
    return Vec2f(v[0] - other.v[0], v[1] - other.v[1]);
} 
Vec2f Vec2f::operator-(double factor) const {
    return Vec2f(v[0] - factor, v[1] - factor);
}
Vec2f Vec2f::operator-() const {
    return Vec2f(-v[0], -v[1]);
} 
bool Vec2f::operator==(const Vec2f &other) const {
    return (v[0] == other.v[0]) && (v[1] == other.v[1]);
}
const Vec2f &Vec2f::operator*=(double scale) {
    v[0] *= scale;
    v[1] *= scale;
    return *this;
} 
const Vec2f &Vec2f::operator*=(const Vec2f &other) {
    v[0] *= other[0];
    v[1] *= other[1];
    return *this;
}
const Vec2f &Vec2f::operator/=(double scale) {
    v[0] /= scale;
    v[1] /= scale;
    return *this;
} 
const Vec2f &Vec2f::operator/=(size_t scale) {
    v[0] /= scale;
    v[1] /= scale;
    return *this;
} 
const Vec2f &Vec2f::operator/=(const Vec2f &other) {
    v[0] /= other[0];
    v[1] /= other[1];
    return *this;
}
const Vec2f &Vec2f::operator+=(const Vec2f &other) {
    v[0] += other.v[0];
    v[1] += other.v[1];
    return *this;
} 
const Vec2f &Vec2f::operator+=(double factor) {
    v[0] += factor;
    v[1] += factor;
    return *this;
} 
const Vec2f &Vec2f::operator-=(const Vec2f &other) {
    v[0] -= other.v[0];
    v[1] -= other.v[1];
    return *this;
} 
const Vec2f &Vec2f::operator-=(double factor) {
    v[0] -= factor;
    v[1] -= factor;
    return *this;
}
Vec2f &Vec2f::min(const Vec2f &other){
    v[0] = std::min(v[0], other[0]);
    v[1] = std::min(v[1], other[1]);
    return *this;
}
Vec2f &Vec2f::min(double other){
    v[0] = std::min(v[0], other);
    v[1] = std::min(v[1], other);
    return *this;
}
Vec2f &Vec2f::max(const Vec2f &other){
    v[0] = std::max(v[0], other[0]);
    v[1] = std::max(v[1], other[1]);
    return *this;
}
Vec2f &Vec2f::max(double other){
    v[0] = std::max(v[0], other);
    v[1] = std::max(v[1], other);
    return *this;
}
Vec2f Vec2f::getMin(const Vec2f &other) const {
    return Vec2f(std::min(v[0], other[0]), std::min(v[1], other[1]));
}
Vec2f Vec2f::getMin(double other) const {
    return Vec2f(std::min(v[0], other), std::min(v[1], other));
}
Vec2f Vec2f::getMax(const Vec2f &other) const{
    return Vec2f(std::max(v[0], other[0]), std::max(v[1], other[1]));
}
Vec2f Vec2f::getMax(double other) const{
    return Vec2f(std::max(v[0], other), std::max(v[1], other));
}
double Vec2f::magnitude() const {
    return std::sqrt(v[0] * v[0] + v[1] * v[1]);
} 
double Vec2f::magnitudeSquared() const {
    return v[0] * v[0] + v[1] * v[1];
} 
Vec2f Vec2f::normalize() const {
    const double m = std::sqrt(v[0] * v[0] + v[1] * v[1]);
    return Vec2f(v[0] / m, v[1] / m);
} 
const Vec2f &Vec2f::normalize_inplace() {
    const double m = std::sqrt(v[0] * v[0] + v[1] * v[1]);
    v[0] /= m;
    v[1] /= m;
    return *this;
}
double Vec2f::dot(const Vec2f &other) const {
    return v[0] * other.v[0] + v[1] * other.v[1];
} 
const Vec2f &Vec2f::to_sph(){ // CHECK outputs nan
    // [r, phi]
    const double temp = std::atan2(v[1], v[0]);
    v[0] = magnitude();
    v[1] = temp;
    return *this;
}
const Vec2f &Vec2f::to_xy(){
    const double temp = v[0];
    v[0] = temp*std::cos(v[1]);
    v[1] = temp*std::sin(v[1]);
    return *this;
}
const Vec2f &Vec2f::to_xy_offset(const Vec2f &ref1, const Vec2f &ref2){
    const Vec2f temp = Vec2f(v[0]*std::cos(v[1]), v[0]*std::sin(v[1])); // CHECK could be better
    v[0] = ref1[0] * temp[0] + ref2[0] * temp[1];
    v[1] = ref1[1] * temp[0] + ref2[1] * temp[1];
    return *this;
}
Vec2f Vec2f::get_sph() const {
    const double r = magnitude();
    return Vec2f(magnitude(), std::atan2(v[1], v[0]));
}
Vec2f Vec2f::get_xy() const {
    return Vec2f(v[0]*std::cos(v[1]), v[0]*std::sin(v[1]));
}
Vec2f Vec2f::get_xy_offset(const Vec2f &ref1, const Vec2f &ref2) const {
    return ref1 * v[0]*std::cos(v[1]) + ref2 * v[0]*std::sin(v[1]);
}
Vec2f Vec2f::ln() const {
    return Vec2f(std::log(v[0]), std::log(v[1]));
}
Vec2f Vec2f::sqrt() const {
    return Vec2f(std::sqrt(v[0]), std::sqrt(v[1]));
}
Vec2f Vec2f::exp() const {
    return Vec2f(std::exp(v[0]), std::exp(v[1]));
}
Vec2f Vec2f::pow(double exp) const {
    return Vec2f(std::pow(v[0], exp), std::pow(v[1], exp));
}
Vec2f &Vec2f::pow_inplace(double exp){
    v[0] = std::pow(v[0], exp);
    v[1] = std::pow(v[1], exp);
    return *this;
}
Vec2f Vec2f::floor() const {
    return Vec2f(std::floor(v[0]), std::floor(v[1]));
}
Vec2f Vec2f::ceil() const {
    return Vec2f(std::ceil(v[0]), std::ceil(v[1]));
}
Vec2f &Vec2f::round_inplace(){
    v[0] = std::round(v[0]);
    v[1] = std::round(v[1]);
    return *this;
}
Vec2f &Vec2f::clamp(double minimum, double maximum){
    min(maximum);
    max(minimum);
    return *this;
}
std::ostream &operator<<(std::ostream &output, const Vec2f &v) {
    std::cout << '[' << v[0] << ", " << v[1] << ']';
    return output;
} 
Vec2f operator*(const double factor, const Vec2f &v) {
    return Vec2f(v[0] * factor, v[1] * factor);
}