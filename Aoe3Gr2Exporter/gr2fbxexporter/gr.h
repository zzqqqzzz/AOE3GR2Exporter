#pragma once
#include <cstdint>
#include <vector>

#include "cgmath.h"
#include "granny.h"
#include "virtual_ptr.h"

enum GR2_curve_format {
    DaKeyframes32f = 0, // Not found in NWN2 files yet
    DaK32fC32f = 1,
    DaIdentity = 2,
    DaConstant32f = 3,
    D3Constant32f = 4,
    D4Constant32f = 5,
    DaK16uC16u = 6,
    DaK8uC8u = 7, // Not found in NWN2 files yet
    D4nK16uC15u = 8,
    D4nK8uC7u = 9,
    D3K16uC16u = 10,
    D3K8uC8u = 11
};

class GR2_DaK16uC16u_view {
public:
    GR2_DaK16uC16u_view(granny_curve_data_da_k16u_c16u& data);

    const std::vector<uint16_t>& encoded_knots() const;
    const std::vector<float>& knots() const;
    const std::vector<uint16_t>& encoded_controls() const;
    const std::vector<float>& controls() const;

private:
    std::vector<uint16_t> encoded_knots_;
    std::vector<float> knots_;
    std::vector<uint16_t> encoded_controls_;
    std::vector<float> controls_;
};

class GR2_DaK32fC32f_view {
public:
    GR2_DaK32fC32f_view(granny_curve_data_da_k32f_c32f& data);

    const std::vector<float>& knots() const;
    const std::vector<Vector4<float>>& controls() const;

private:
    std::vector<float> knots_;
    std::vector<Vector4<float>> controls_;
};

class GR2_D4nK16uC15u_view {
public:
    uint16_t selectors[4];
    float scales[4];
    float offsets[4];

    GR2_D4nK16uC15u_view(granny_curve_data_d4n_k16u_c15u& data);

    const std::vector<uint16_t>& encoded_knots() const;
    const std::vector<float>& knots() const;
    const std::vector<Vector3<uint16_t>>& encoded_controls() const;
    const std::vector<Vector4<float>>& controls() const;

private:
    std::vector<uint16_t> encoded_knots_;
    std::vector<float> knots_;
    std::vector<Vector3<uint16_t>> encoded_controls_;
    std::vector<Vector4<float>> controls_;
};


class GR2_D4nK8uC7u_view {
public:
    uint16_t selectors[4];
    float scales[4];
    float offsets[4];

    GR2_D4nK8uC7u_view(granny_curve_data_d4n_k8u_c7u& data);

    const std::vector<uint8_t>& encoded_knots() const;
    const std::vector<float>& knots() const;
    const std::vector<Vector3<uint8_t>>& encoded_controls() const;
    const std::vector<Vector4<float>>& controls() const;

private:
    std::vector<uint8_t> encoded_knots_;
    std::vector<float> knots_;
    std::vector<Vector3<uint8_t>> encoded_controls_;
    std::vector<Vector4<float>> controls_;
};


class GR2_D3K16uC16u_view {
public:
    GR2_D3K16uC16u_view(granny_curve_data_d3_k16u_c16u& data);

    const std::vector<uint16_t>& encoded_knots() const;
    const std::vector<float>& knots() const;
    const std::vector<Vector3<uint16_t>>& encoded_controls() const;
    const std::vector<Vector3<float>>& controls() const;

private:
    std::vector<uint16_t> encoded_knots_;
    std::vector<float> knots_;
    std::vector<Vector3<uint16_t>> encoded_controls_;
    std::vector<Vector3<float>> controls_;
};


class GR2_D3K8uC8u_view {
public:
    GR2_D3K8uC8u_view(granny_curve_data_d3_k8u_c8u& data);

    const std::vector<uint8_t>& encoded_knots() const;
    const std::vector<float>& knots() const;
    const std::vector<Vector3<uint8_t>>& encoded_controls() const;
    const std::vector<Vector3<float>>& controls() const;

private:
    std::vector<uint8_t> encoded_knots_;
    std::vector<float> knots_;
    std::vector<Vector3<uint8_t>> encoded_controls_;
    std::vector<Vector3<float>> controls_;
};

class GR2_curve_view {
public:
    GR2_curve_view(granny_curve2& curve);

    uint8_t degree() const;
    const std::vector<float>& knots() const;
    const std::vector<Vector4<float>>& controls() const;

private:
    uint8_t degree_;
    std::vector<float> knots_;
    std::vector<Vector4<float>> controls_;
};