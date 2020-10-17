#include <cmath>
#include <string.h>
#include "gr.h"

static const float scale_table[] = {
    1.4142135f,   0.70710677f,  0.35355338f,  0.35355338f,
    0.35355338f,  0.17677669f,  0.17677669f,  0.17677669f,
    -1.4142135f,  -0.70710677f, -0.35355338f, -0.35355338f,
    -0.35355338f, -0.17677669f, -0.17677669f, -0.17677669f };

static const float offset_table[] = {
    -0.70710677f, -0.35355338f, -0.53033006f,  -0.17677669f,
    0.17677669f,  -0.17677669f, -0.088388346f, 0.0f,
    0.70710677f,  0.35355338f,  0.53033006f,   0.17677669f,
    -0.17677669f, 0.17677669f,  0.088388346f,  -0.0f };

using namespace std;

static void compute_offsets(float offsets[4], uint16_t selectors[4])
{
    for (int i = 0; i < 4; ++i)
        offsets[i] = offset_table[selectors[i]];
}

static void compute_selectors(uint16_t selectors[4],
    uint16_t scale_offset_table_entries)
{
    selectors[0] = (scale_offset_table_entries >> 0) & 0x0F;
    selectors[1] = (scale_offset_table_entries >> 4) & 0x0F;
    selectors[2] = (scale_offset_table_entries >> 8) & 0x0F;
    selectors[3] = (scale_offset_table_entries >> 12) & 0x0F;
}


static Vector4<float> decode_D4nK16uC15u(uint16_t a, uint16_t b, uint16_t c,
    float scales[], float offsets[])
{
    // A quaternion (4 components) is encoded in three values (a, b, c)
    //
    // a: 15 ... 1 0 | b: 15 ... 1 0 | c: 15 ... 1 0
    //    g    da        s1a   db        s1b   dc
    //
    // da, db, dc: 3 components of the quaternion
    // g: sign flag for 4th component of the quaternion (dd)
    // s1a, s1b: swizzle

    int s1a = (b & 0x8000) >> 14;
    int s1b = c >> 15;
    int swizzle1 = s1a | s1b;

    // swizzle_n = swizzle_{n-1} mod 4
    int swizzle2 = (swizzle1 + 1) & 3;
    int swizzle3 = (swizzle2 + 1) & 3;
    int swizzle4 = (swizzle3 + 1) & 3;

    float da = (a & 0x7fff) * scales[swizzle2] + offsets[swizzle2];
    float db = (b & 0x7fff) * scales[swizzle3] + offsets[swizzle3];
    float dc = (c & 0x7fff) * scales[swizzle4] + offsets[swizzle4];

    // Reconstruct dd considering quaternion is unit length
    float dd = sqrtf(1 - (da * da + db * db + dc * dc));
    if ((a & 0x8000) != 0)
        dd = -dd;

    Vector4<float> quat;
    quat[swizzle2] = da;
    quat[swizzle3] = db;
    quat[swizzle4] = dc;
    quat[swizzle1] = dd;
    return quat;
}

static Vector4<float> decode_D4nK8uC7u(uint8_t a, uint8_t b, uint8_t c,
    float scales[], float offsets[])
{
    // A quaternion (4 components) is encoded in three values (a, b, c)
    //
    // a: 7 ... 1 0 | b: 7 ... 1 0 | c: 7 ... 1 0
    //    g   da       s1a   db       s1b   dc
    //
    // da, db, dc: 3 of 4 components
    // g: sign flag for 4th component (dd)
    // s1a, s1b: swizzle

    int s1a = (b & 0x80) >> 6;
    int s1b = (c & 0x80) >> 7;
    int swizzle1 = s1a | s1b;

    // swizzle_n = swizzle_{n-1} mod 4
    int swizzle2 = (swizzle1 + 1) & 3;
    int swizzle3 = (swizzle2 + 1) & 3;
    int swizzle4 = (swizzle3 + 1) & 3;

    float da = (a & 0x7f) * scales[swizzle2] + offsets[swizzle2];
    float db = (b & 0x7f) * scales[swizzle3] + offsets[swizzle3];
    float dc = (c & 0x7f) * scales[swizzle4] + offsets[swizzle4];

    // Reconstruct dd considering quaternion is unit length
    float dd = sqrtf(1 - (da * da + db * db + dc * dc));
    if ((a & 0x80) != 0)
        dd = -dd;

    Vector4<float> quat;
    quat[swizzle2] = da;
    quat[swizzle3] = db;
    quat[swizzle4] = dc;
    quat[swizzle1] = dd;
    return quat;
}

GR2_DaK16uC16u_view::GR2_DaK16uC16u_view(granny_curve_data_da_k16u_c16u& data)
{
    int dimension = data.ControlScaleOffsetCount / 2;
    int knots_count = data.KnotControlCount / (dimension + 1);

    float one_over_knot_scale;
    unsigned tmp = (unsigned)data.OneOverKnotScaleTrunc << 16;
    memcpy(&one_over_knot_scale, &tmp, sizeof(tmp));

    for (int i = 0; i < knots_count; ++i) {
        encoded_knots_.push_back(data.KnotsControls[i]);
        knots_.push_back(data.KnotsControls[i] / one_over_knot_scale);
    }

    int controls_count = data.KnotControlCount - knots_count;
    auto controls = data.KnotsControls + knots_count;

    for (int i = 0; i < controls_count; ++i) {
        encoded_controls_.push_back(controls[i]);
        int scale_index = i % dimension;
        int offset_index = scale_index + dimension;
        controls_.push_back(controls[i] * data.ControlScaleOffsets[scale_index] + data.ControlScaleOffsets[offset_index]);
    }
}

const std::vector<uint16_t>& GR2_DaK16uC16u_view::encoded_knots() const
{
    return encoded_knots_;
}

const std::vector<float>& GR2_DaK16uC16u_view::knots() const
{
    return knots_;
}

const std::vector<uint16_t>& GR2_DaK16uC16u_view::encoded_controls() const
{
    return encoded_controls_;
}

const std::vector<float>& GR2_DaK16uC16u_view::controls() const
{
    return controls_;
}

GR2_DaK32fC32f_view::GR2_DaK32fC32f_view(granny_curve_data_da_k32f_c32f& data)
{
    for (int i = 0; i < data.KnotCount; ++i)
        knots_.push_back(data.Knots[i]);

    if (data.KnotCount * 3 == data.ControlCount) { // Positions
        for (int i = 0; i < data.KnotCount; ++i) {
            float x = data.Controls[i * 3 + 0];
            float y = data.Controls[i * 3 + 1];
            float z = data.Controls[i * 3 + 2];
            controls_.emplace_back(x, y, z, 1.0f);
        }
    }
    else if (data.KnotCount * 4 == data.ControlCount) { // Quaternions
        for (int i = 0; i < data.KnotCount; ++i) {
            float x = data.Controls[i * 4 + 0];
            float y = data.Controls[i * 4 + 1];
            float z = data.Controls[i * 4 + 2];
            float w = data.Controls[i * 4 + 3];
            controls_.emplace_back(x, y, z, w);
        }
    }
}

const std::vector<float>& GR2_DaK32fC32f_view::knots() const
{
    return knots_;
}

const std::vector<Vector4<float>>& GR2_DaK32fC32f_view::controls() const
{
    return controls_;
}

GR2_D4nK16uC15u_view::GR2_D4nK16uC15u_view(granny_curve_data_d4n_k16u_c15u& data)
{
    int knots_count = data.KnotControlCount / 4;
    for (int i = 0; i < knots_count; ++i) {
        encoded_knots_.push_back(data.KnotsControls[i]);
        knots_.push_back(data.KnotsControls[i] /
            data.OneOverKnotScale);
    }

    compute_selectors(selectors, data.ScaleOffsetTableEntries);

    for (int i = 0; i < 4; ++i)
        scales[i] = scale_table[selectors[i]] * 0.000030518509f;

    compute_offsets(offsets, selectors);

    uint16_t * controls = data.KnotsControls + knots_count;
    int controls_count = data.KnotControlCount - knots_count;

    for (int i = 0; i < controls_count; i += 3) {
        uint16_t a = controls[i + 0];
        uint16_t b = controls[i + 1];
        uint16_t c = controls[i + 2];
        encoded_controls_.emplace_back(a, b, c);
        controls_.push_back(
            decode_D4nK16uC15u(a, b, c, scales, offsets));
    }
}

const std::vector<uint16_t>& GR2_D4nK16uC15u_view::encoded_knots() const
{
    return encoded_knots_;
}

const std::vector<float>& GR2_D4nK16uC15u_view::knots() const
{
    return knots_;
}

const std::vector<Vector3<uint16_t>>&
GR2_D4nK16uC15u_view::encoded_controls() const
{
    return encoded_controls_;
}

const std::vector<Vector4<float>>& GR2_D4nK16uC15u_view::controls() const
{
    return controls_;
}


GR2_D4nK8uC7u_view::GR2_D4nK8uC7u_view(granny_curve_data_d4n_k8u_c7u& data)
{
    int knots_count = data.KnotControlCount / 4;
    for (int i = 0; i < knots_count; ++i) {
        encoded_knots_.push_back(data.KnotsControls[i]);
        knots_.push_back(data.KnotsControls[i] /
            data.OneOverKnotScale);
    }

    compute_selectors(selectors, data.ScaleOffsetTableEntries);

    for (int i = 0; i < 4; ++i)
        scales[i] = scale_table[selectors[i]] * 0.0078740157f;

    compute_offsets(offsets, selectors);

    uint8_t * controls = data.KnotsControls + knots_count;
    int controls_count = data.KnotControlCount - knots_count;

    for (int i = 0; i < controls_count; i += 3) {
        uint8_t a = controls[i + 0];
        uint8_t b = controls[i + 1];
        uint8_t c = controls[i + 2];
        encoded_controls_.emplace_back(a, b, c);
        controls_.emplace_back(
            decode_D4nK8uC7u(a, b, c, scales, offsets));
    }
}

const std::vector<uint8_t>& GR2_D4nK8uC7u_view::encoded_knots() const
{
    return encoded_knots_;
}

const std::vector<float>& GR2_D4nK8uC7u_view::knots() const
{
    return knots_;
}

const std::vector<Vector3<uint8_t>>&
GR2_D4nK8uC7u_view::encoded_controls() const
{
    return encoded_controls_;
}

const std::vector<Vector4<float>>& GR2_D4nK8uC7u_view::controls() const
{
    return controls_;
}


GR2_D3K16uC16u_view::GR2_D3K16uC16u_view(granny_curve_data_d3_k16u_c16u& data)
{
    int knots_count = data.KnotControlCount / 4;

    float one_over_knot_scale;
    unsigned tmp = (unsigned)data.OneOverKnotScaleTrunc << 16;
    memcpy(&one_over_knot_scale, &tmp, sizeof(tmp));

    for (int i = 0; i < knots_count; ++i) {
        encoded_knots_.push_back(data.KnotsControls[i]);
        knots_.push_back(data.KnotsControls[i] / one_over_knot_scale);
    }

    int controls_count = data.KnotControlCount - knots_count;
    auto controls = data.KnotsControls + knots_count;

    for (int i = 0; i < controls_count; i += 3) {
        auto a = controls[i + 0];
        auto b = controls[i + 1];
        auto c = controls[i + 2];
        encoded_controls_.emplace_back(a, b, c);

        float x = a * data.ControlScales[0] + data.ControlScales[0];
        float y = b * data.ControlScales[1] + data.ControlScales[1];
        float z = c * data.ControlScales[2] + data.ControlScales[2];
        controls_.emplace_back(x, y, z);
    }
}

const std::vector<uint16_t>& GR2_D3K16uC16u_view::encoded_knots() const
{
    return encoded_knots_;
}

const std::vector<float>& GR2_D3K16uC16u_view::knots() const
{
    return knots_;
}

const std::vector<Vector3<uint16_t>>& GR2_D3K16uC16u_view::encoded_controls() const
{
    return encoded_controls_;
}

const std::vector<Vector3<float>>& GR2_D3K16uC16u_view::controls() const
{
    return controls_;
}

GR2_D3K8uC8u_view::GR2_D3K8uC8u_view(granny_curve_data_d3_k8u_c8u& data)
{
    int knots_count = data.KnotControlCount / 4;

    float one_over_knot_scale;
    unsigned tmp = (unsigned)data.OneOverKnotScaleTrunc << 16;
    memcpy(&one_over_knot_scale, &tmp, sizeof(tmp));

    for (int i = 0; i < knots_count; ++i) {
        encoded_knots_.push_back(data.KnotsControls[i]);
        knots_.push_back(data.KnotsControls[i] / one_over_knot_scale);
    }

    int controls_count = data.KnotControlCount - knots_count;
    uint8_t* controls = data.KnotsControls + knots_count;

    for (int i = 0; i < controls_count; i += 3) {
        uint8_t a = controls[i + 0];
        uint8_t b = controls[i + 1];
        uint8_t c = controls[i + 2];
        encoded_controls_.emplace_back(a, b, c);

        float x = a * data.ControlScales[0] + data.ControlScales[0];
        float y = b * data.ControlScales[1] + data.ControlScales[1];
        float z = c * data.ControlScales[2] + data.ControlScales[2];
        controls_.emplace_back(x, y, z);
    }
}

const std::vector<uint8_t>& GR2_D3K8uC8u_view::encoded_knots() const
{
    return encoded_knots_;
}

const std::vector<float>& GR2_D3K8uC8u_view::knots() const
{
    return knots_;
}

const std::vector<Vector3<uint8_t>>& GR2_D3K8uC8u_view::encoded_controls() const
{
    return encoded_controls_;
}

const std::vector<Vector3<float>>& GR2_D3K8uC8u_view::controls() const
{
    return controls_;
}

GR2_curve_view::GR2_curve_view(granny_curve2& curve)
{
    granny_curve_data_header* header = (granny_curve_data_header*)curve.CurveData.Object;
    degree_ = header->Degree;

    if (header->Format == DaK32fC32f) {
        granny_curve_data_da_k32f_c32f* data =
            (granny_curve_data_da_k32f_c32f*)curve.CurveData.Object;
        GR2_DaK32fC32f_view view(*data);
        knots_ = view.knots();
        controls_ = view.controls();
    }
    else if (header->Format == DaIdentity) {
    }
    else if (header->Format == D3Constant32f) {
        granny_curve_data_da_constant32f* data =
            (granny_curve_data_da_constant32f*)curve.CurveData.Object;
        knots_.push_back(0.0f);
        controls_.emplace_back(data->Controls[0],
            data->Controls[1], data->Controls[2], 1.0f);
    }
    else if (header->Format == D4nK16uC15u) {
        granny_curve_data_d4n_k16u_c15u* data =
            (granny_curve_data_d4n_k16u_c15u*)curve.CurveData.Object;
        GR2_D4nK16uC15u_view view(*data);
        knots_ = view.knots();
        controls_ = view.controls();
    }
    else if (header->Format == D4nK8uC7u) {
        granny_curve_data_d4n_k8u_c7u* data =
            (granny_curve_data_d4n_k8u_c7u*)curve.CurveData.Object;
        GR2_D4nK8uC7u_view view(*data);
        knots_ = view.knots();
        controls_ = view.controls();
    }
    else if (header->Format == D3K16uC16u) {
        granny_curve_data_d3_k16u_c16u* data =
            (granny_curve_data_d3_k16u_c16u*)curve.CurveData.Object;
        GR2_D3K16uC16u_view view(*data);
        knots_ = view.knots();
        for (auto& c : view.controls())
            controls_.emplace_back(c.x, c.y, c.z, 1.0f);
    }
    else if (header->Format == D3K8uC8u) {
        granny_curve_data_d3_k8u_c8u* data =
            (granny_curve_data_d3_k8u_c8u*)curve.CurveData.Object;
        GR2_D3K8uC8u_view view(*data);
        knots_ = view.knots();
        for (auto& c : view.controls())
            controls_.emplace_back(c.x, c.y, c.z, 1.0f);
    }
}

uint8_t GR2_curve_view::degree() const
{
    return degree_;
}

const std::vector<float>& GR2_curve_view::knots() const
{
    return knots_;
}

const std::vector<Vector4<float>>& GR2_curve_view::controls() const
{
    return controls_;
}
