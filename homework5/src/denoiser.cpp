#include "denoiser.h"

Denoiser::Denoiser() : m_useTemportal(false) {}

void Denoiser::Reprojection(const FrameInfo &frameInfo) {
    int height = m_accColor.m_height;
    int width = m_accColor.m_width;
    Matrix4x4 preWorldToScreen =
        m_preFrameInfo.m_matrix[m_preFrameInfo.m_matrix.size() - 1];
    Matrix4x4 preWorldToCamera =
        m_preFrameInfo.m_matrix[m_preFrameInfo.m_matrix.size() - 2];
#pragma omp parallel for
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // TODO: Reproject
            auto ipos = frameInfo.m_position(x, y);

            auto iobjid = frameInfo.m_id(x, y);
            if (iobjid < 0) {
                m_valid(x, y) = false;
                m_misc(x, y) = Float3(0.f);
                continue;
            }
            auto iobj2world = frameInfo.m_matrix[iobjid];
            auto lastObj2World = m_preFrameInfo.m_matrix[iobjid];
            auto lastUV = preWorldToScreen(
                lastObj2World(Inverse(iobj2world)(ipos, Float3::EType::Point),
                              Float3::EType::Point),
                Float3::EType::Point);
            if (lastUV.x < 0 || lastUV.x > width || lastUV.y < 0 || lastUV.y > height) {
                m_valid(x, y) = false;
                m_misc(x, y) = Float3(0.f);
            } else {
                int lastX = lastUV.x;
                int lastY = lastUV.y;

                auto lastUVObjId = m_preFrameInfo.m_id(lastX, lastY);
                if (lastUVObjId != iobjid) {
                    m_valid(x, y) = false;
                    m_misc(x, y) = Float3(0.f);
                } else {
                    m_valid(x, y) = true;
                    m_misc(x, y) = m_accColor(lastX, lastY);//here use had filtered cache color
                }                
            }
            //m_valid(x, y) = false;
            //m_misc(x, y) = Float3(0.f);
        }
    }
    std::swap(m_misc, m_accColor);
}

Float3 Denoiser::ClampAcc(int x, int y, int kernelRadius) {
    int height = m_accColor.m_height;
    int width = m_accColor.m_width;

    int xmin = std::max(0, x - kernelRadius);
    int xmax = std::min(width, kernelRadius + x);

    int ymin = std::max(0, y - kernelRadius);
    int ymax = std::min(height, kernelRadius + y);

    Float3 e;
    Float3 u; 
    Float3 color = m_accColor(x, y);

    int validCount = 0;
    for (int i = xmin; i < xmax; i++) {
        for (int j = ymin; j < ymax; j++) {
            if (m_valid(i, j)) {
                Float3 color = m_accColor(i, j);
                e += color;
                validCount++;
            }
        }
    }

    e = e / validCount;
    for (int i = xmin; i < xmax; i++) {
        for (int j = ymin; j < ymax; j++) {
            if (m_valid(i, j)) {
                Float3 color = m_accColor(i, j);
                u += (e - color) * (e - color);
            }
        }
    }
    u = u / validCount;
    color = Clamp(color, e - u * m_colorBoxK, e + u * m_colorBoxK);
    return color;
}

void Denoiser::TemporalAccumulation(const Buffer2D<Float3> &curFilteredColor) {
    int height = m_accColor.m_height;
    int width = m_accColor.m_width;
    int kernelRadius = 3;
#pragma omp parallel for
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // TODO: Temporal clamp
            Float3 color = m_accColor(x, y);
            //Float3 color = ClampAcc(x, y, kernelRadius);


            // TODO: Exponential moving average
            float alpha = 1.0;
            if (m_valid(x, y)) {
                alpha = m_alpha;
                color = ClampAcc(x, y, kernelRadius);
            }            
            m_misc(x, y) = Lerp(color, curFilteredColor(x, y), alpha);
        }
    }
    std::swap(m_misc, m_accColor);
}


Float3 Denoiser::FilterPixel(const FrameInfo& frameInfo, int x, int y, int kernelRadius) {
    int height = frameInfo.m_beauty.m_height;
    int width = frameInfo.m_beauty.m_width;
    int xmin = std::max(0, x - kernelRadius);
    int xmax = std::min(width, kernelRadius + x);

    int ymin = std::max(0, y - kernelRadius);
    int ymax = std::min(height, kernelRadius + y);

    Float3 split_sum_weigth;
    Float3 split_sum;
    Float3 iColor = frameInfo.m_beauty(x, y);
    Float3 iPos = frameInfo.m_position(x, y);
    Float3 iNormal = frameInfo.m_normal(x, y);
    if (iNormal.x == 0 && iNormal.y == 0 && iNormal.z == 0) {
        return frameInfo.m_beauty(x, y);
    }
    //FLOAT iObjId = frameInfo.m_id(x, y);
    //float* weights = new float[xmax - xmin, ymax - ymin];
    //Float3 *splitColors = new Float3[xmax - xmin, ymax - ymin];
    //if (x == 564 && y == 300) {
    //    int dsadasd = 0;
    //}
    for (int i = xmin; i < xmax; i++) {
        for (int j = ymin; j < ymax; j++) {
            if (i == x && j == y) {
                split_sum_weigth += 1;
                split_sum += frameInfo.m_beauty(i, j) * 1;
            } else {
                Float3 jColor = frameInfo.m_beauty(i, j);
                Float3 jPos = frameInfo.m_position(i, j);
                Float3 jNormal = frameInfo.m_normal(i, j);
                if (jNormal.x == 0 && jNormal.y == 0 && jNormal.z == 0) {
                    continue;
                }
                //FLOAT jObjId = frameInfo.m_id(i, j);
                float dNormal = SafeAcos(Dot(iNormal, jNormal));
                auto dir = jPos - iPos;
                float dPlane = Dot(iNormal, Normalize(dir));
                float dDis = (x - i) * (x - i) + (y - j) * (y - j);
                auto iL = Luminance(iColor);
                auto jL = Luminance(jColor);
                float dCor = abs( iL - jL);

                Float3 weigth =
                    exp(-dDis / 2.0 / m_sigmaCoord / m_sigmaCoord -
                        dCor * dCor / 2.0 / m_sigmaColor / m_sigmaColor -
                        dNormal * dNormal / 2.0 / m_sigmaNormal / m_sigmaNormal -
                        dPlane * dPlane / 2.0 / m_sigmaPlane / m_sigmaPlane);
                Float3 weigthColor = frameInfo.m_beauty(i, j) * weigth;                
                split_sum_weigth += weigth;
                split_sum += weigthColor;
                //weights[(i - xmin) , (j - ymin)] = weigth.x;
                //splitColors[(i - xmin), (j - ymin)] = weigthColor;
                //if (x == 564 && y == 300) {
                //    int dsadasd = 0;
                //}
            }

        }
    }
    //if (x == 564 && y == 300) {
    //    int dsadasd = 0;
    //}
    auto ret = split_sum / split_sum_weigth;
    return ret;
}

Buffer2D<Float3> Denoiser::Filter(const FrameInfo &frameInfo) {
    int height = frameInfo.m_beauty.m_height;
    int width = frameInfo.m_beauty.m_width;
    Buffer2D<Float3> filteredImage = CreateBuffer2D<Float3>(width, height);
    int kernelRadius = 16;
#pragma omp parallel for
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // TODO: Joint bilateral filter                        
            filteredImage(x, y) = FilterPixel(frameInfo,x,y,kernelRadius);
            //filteredImage(x, y) = frameInfo.m_beauty(x, y);
        }
    }
    return filteredImage;
}

void Denoiser::Init(const FrameInfo &frameInfo, const Buffer2D<Float3> &filteredColor) {
    m_accColor.Copy(filteredColor);
    int height = m_accColor.m_height;
    int width = m_accColor.m_width;
    m_misc = CreateBuffer2D<Float3>(width, height);
    m_valid = CreateBuffer2D<bool>(width, height);
}

void Denoiser::Maintain(const FrameInfo &frameInfo) { m_preFrameInfo = frameInfo; }

Buffer2D<Float3> Denoiser::ProcessFrame(const FrameInfo &frameInfo) {
    // Filter current frame
    Buffer2D<Float3> filteredColor;
    filteredColor = Filter(frameInfo);

    // Reproject previous frame color to current
    if (m_useTemportal) {
        Reprojection(frameInfo);
        TemporalAccumulation(filteredColor);
    } else {
        Init(frameInfo, filteredColor);
    }

    // Maintain
    Maintain(frameInfo);
    if (!m_useTemportal) {
        m_useTemportal = true;
    }
    return m_accColor;
}
