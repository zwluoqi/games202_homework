#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/ray.h>
#include <filesystem/resolver.h>
#include <sh/spherical_harmonics.h>
#include <sh/default_image.h>
#include <Eigen/Core>
#include <fstream>
#include <random>
#include <stb_image.h>

NORI_NAMESPACE_BEGIN

namespace ProjEnv
{
    std::vector<std::unique_ptr<float[]>>
    LoadCubemapImages(const std::string &cubemapDir, int &width, int &height,
                      int &channel)
    {
        std::vector<std::string> cubemapNames{"negx.jpg", "posx.jpg", "posy.jpg",
                                              "negy.jpg", "posz.jpg", "negz.jpg"};
        std::vector<std::unique_ptr<float[]>> images(6);
        for (int i = 0; i < 6; i++)
        {
            std::string filename = cubemapDir + "/" + cubemapNames[i];
            int w, h, c;
            float *image = stbi_loadf(filename.c_str(), &w, &h, &c, 3);
            if (!image)
            {
                std::cout << "Failed to load image: " << filename << std::endl;
                exit(-1);
            }
            if (i == 0)
            {
                width = w;
                height = h;
                channel = c;
            }
            else if (w != width || h != height || c != channel)
            {
                std::cout << "Dismatch resolution for 6 images in cubemap" << std::endl;
                exit(-1);
            }
            images[i] = std::unique_ptr<float[]>(image);
            int index = (0 * 128 + 0) * channel;
            // std::cout << images[i][index + 0] << "\t" << images[i][index + 1] << "\t"
            //           << images[i][index + 2] << std::endl;
        }
        return images;
    }

    const Eigen::Vector3f cubemapFaceDirections[6][3] = {
        {{0, 0, 1}, {0, -1, 0}, {-1, 0, 0}},  // negx
        {{0, 0, 1}, {0, -1, 0}, {1, 0, 0}},   // posx
        {{1, 0, 0}, {0, 0, -1}, {0, -1, 0}},  // negy
        {{1, 0, 0}, {0, 0, 1}, {0, 1, 0}},    // posy
        {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}, // negz
        {{1, 0, 0}, {0, -1, 0}, {0, 0, 1}},   // posz
    };

    struct SHSample {
        Eigen::Vector3d sph;
        Eigen::Vector3d vec;
        Eigen::Array3f coeff;
    };

    struct InterRef {
        Point3f p;
        Normal3f n;
        int vertexIndex;
    };

    float CalcPreArea(const float &x, const float &y)
    {
        return std::atan2(x * y, std::sqrt(x * x + y * y + 1.0));
    }

    float CalcArea(const float &u_, const float &v_, const int &width,
                   const int &height)
    {
        // transform from [0..res - 1] to [- (1 - 1 / res) .. (1 - 1 / res)]
        // ( 0.5 is for texel center addressing)
        float u = (2.0 * (u_ + 0.5) / width) - 1.0;
        float v = (2.0 * (v_ + 0.5) / height) - 1.0;

        // shift from a demi texel, mean 1.0 / size  with u and v in [-1..1]
        float invResolutionW = 1.0 / width;
        float invResolutionH = 1.0 / height;

        // u and v are the -1..1 texture coordinate on the current face.
        // get projected area for this texel
        float x0 = u - invResolutionW;
        float y0 = v - invResolutionH;
        float x1 = u + invResolutionW;
        float y1 = v + invResolutionH;
        float angle = CalcPreArea(x0, y0) - CalcPreArea(x0, y1) -
                      CalcPreArea(x1, y0) + CalcPreArea(x1, y1);

        return angle;
    }

    // template <typename T> T ProjectSH() {}

    template <size_t SHOrder>
    std::vector<Eigen::Array3f> PrecomputeCubemapSH(const std::vector<std::unique_ptr<float[]>> &images,
                                                    const int &width, const int &height,
                                                    const int &channel)
    {
        double PI = 3.1415926;

        std::vector<Eigen::Vector3f> cubemapDirs;
       
        cubemapDirs.reserve(6 * width * height);
        for (int i = 0; i < 6; i++)
        {
            Eigen::Vector3f faceDirX = cubemapFaceDirections[i][0];
            Eigen::Vector3f faceDirY = cubemapFaceDirections[i][1];
            Eigen::Vector3f faceDirZ = cubemapFaceDirections[i][2];
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    float u = 2 * ((x + 0.5) / width) - 1;
                    float v = 2 * ((y + 0.5) / height) - 1;
                    Eigen::Vector3f dir = (faceDirX * u + faceDirY * v + faceDirZ).normalized();
                    cubemapDirs.push_back(dir);
                }
            }
        }
        constexpr int SHNum = (SHOrder + 1) * (SHOrder + 1);
        std::vector<Eigen::Array3f> SHCoeffiecents(SHNum);
        for (int i = 0; i < SHNum; i++)
            SHCoeffiecents[i] = Eigen::Array3f(0);
        float sumWeiht = 0;
        for (int i = 0; i < 6; i++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    // TODO: here you need to compute light sh of each face of cubemap of each pixel
                    // TODO: 此处你需要计算每个像素下cubemap某个面的球谐系数
                    Eigen::Vector3f dir = cubemapDirs[i * width * height + y * width + x];
                    int index = (y * width + x) * channel;
                    Eigen::Array3f Le(images[i][index + 0], images[i][index + 1],
                                      images[i][index + 2]);   

                    
                    
                    float gg = CalcArea(x, y, width, height);

                    for (int l = 0; l < SHOrder+1; ++l) {
                        for (int m = -l; m <= l; ++m) {
                            int index = l * (l + 1) + m;
                            const Eigen::Vector3d dirNor = Eigen::Vector3d(dir[0],dir[1],dir[2]);
                            double ci = sh::EvalSH(l, m, dirNor.normalized());
                            SHCoeffiecents[index] += Le * (gg * ci);
                        }
                    }
                    
                }
            }
        }
        //const double area = 4.0 * PI;
        //double factor = area / (6*width*height);
        //for (int i = 0; i < SHNum; ++i) {
        //    SHCoeffiecents[i] = SHCoeffiecents[i] * factor; // NOTE: vector-scalar multiply
        //}
        return SHCoeffiecents;
    }
}

class PRTIntegrator : public Integrator
{
public:
    static constexpr int SHOrder = 2;
    static constexpr int SHCoeffLength = (SHOrder + 1) * (SHOrder + 1);

    enum class Type
    {
        Unshadowed = 0,
        Shadowed = 1,
        Interreflection = 2
    };

    PRTIntegrator(const PropertyList &props)
    {
        /* No parameters this time */
        m_SampleCount = props.getInteger("PRTSampleCount", 100);
        m_CubemapPath = props.getString("cubemap");
        auto type = props.getString("type", "unshadowed");
        if (type == "unshadowed")
        {
            m_Type = Type::Unshadowed;
        }
        else if (type == "shadowed")
        {
            m_Type = Type::Shadowed;
        }
        else if (type == "interreflection")
        {
            m_Type = Type::Interreflection;
            m_Bounce = props.getInteger("bounce", 1);
        }
        else
        {
            throw NoriException("Unsupported type: %s.", type);
        }
    }



    void  interreflectionFunc(int vertexIndex,const Scene* scene, const Point3f& v, const Normal3f& n, int bounces,int sample_count,int order,std::vector<Eigen::MatrixXf> &sh_buffer, std::vector<ProjEnv::InterRef>&  hit) {

        const auto mesh = scene->getMeshes()[0];
        const MatrixXf& normals = mesh->getVertexNormals();


        // This is the approach demonstrated in [1] and is useful for arbitrary
        // functions on the sphere that are represented analytically.
        const int sample_side = static_cast<int>(floor(sqrt(sample_count)));
        std::unique_ptr<std::vector<double>> coeffs(new std::vector<double>());
        coeffs->assign(GetCoefficientCount(order), 0.0);

        // generate sample_side^2 uniformly and stratified samples over the sphere
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> rng(0.0, 1.0);

        for (int t = 0; t < sample_side; t++) {
            for (int p = 0; p < sample_side; p++) {
                double alpha = (t + rng(gen)) / sample_side;
                double beta = (p + rng(gen)) / sample_side;
                // See http://www.bogotobogo.com/Algorithms/uniform_distribution_sphere.php
                double phi = 2.0 * M_PI * beta;
                double theta = acos(2.0 * alpha - 1.0);

                // evaluate the analytic function for the current spherical coords
                //double func_value = func(phi, theta);
                Eigen::Array3d d = sh::ToVector(phi, theta);
                const Vector3f& wi = Vector3f(d.x(), d.y(), d.z()).normalized();
                double s = 0;
                float m = 0.0f;
                m = n.dot(wi);
                double md = (double)m;
                s = std::max(0.0, md);
                if (md > 0.0) {
                    Ray3f ray = Ray3f(v, wi);
                    Intersection its;
                    bool cross = scene->rayIntersect(ray, its);
                    if (cross) {
                        const Eigen::Matrix<Vector3f::Scalar, SHCoeffLength, 1> sh0 = sh_buffer[bounces-1].col(its.tri_index.x()),
                            sh1 = sh_buffer[bounces - 1].col(its.tri_index.y()),
                            sh2 = sh_buffer[bounces - 1].col(its.tri_index.z());
                        const Vector3f& bary = its.bary;
                        const Eigen::Matrix<Vector3f::Scalar, SHCoeffLength, 1> barySH = bary.x() * sh0 + bary.y() * sh1 + bary.z() * sh2;
                        for (int j = 0; j < SHCoeffLength; j++)
                        {
                            (*coeffs)[j] += barySH.coeffRef(j)* md;;
                        }

                        if (bounces < m_Bounce) {
                            const Point3f& p = its.p;
                            const Normal3f& n1 = normals.col(its.tri_index.x()),
                                n2 = normals.col(its.tri_index.y()),
                                n3 = normals.col(its.tri_index.z());
                            const Normal3f& normal = (bary.x() * n1 + bary.y() * n2 + bary.z() * n3).normalized();                           
                            //interreflectionFunc(vertexIndex, scene, p, normal, bounces - 1, sample_count, order, sh_buffer);
                            ProjEnv::InterRef interRef;
                            interRef.p = p;
                            interRef.n = normal;
                            interRef.vertexIndex = vertexIndex;
                            hit.push_back(interRef);
                            //hit_points.push_back(p);
                            //hit_normals.push_back(normal);
                        }
                    }
                }
            }
        }

        // scale by the probability of a particular sample, which is
        // 4pi/sample_side^2. 4pi for the surface area of a unit sphere, and
        // 1/sample_side^2 for the number of samples drawn uniformly.
        double weight = 4.0 * M_PI / (sample_side * sample_side);
        for (int i = 0; i < bounces; i++) {
            weight *= (1/ (4.0 * M_PI));
        }
        for (unsigned int i = 0; i < coeffs->size(); i++) {
            (*coeffs)[i] *= weight;
        }

        for (int j = 0; j < SHCoeffLength; j++)
        {
            sh_buffer[bounces].col(vertexIndex).coeffRef(j) = (*coeffs)[j];
        }
        return ;
        //if (bounces < m_Bounce) {
        //    for (int i = 0; i < hit_points.size(); i++) {
        //        interreflectionFunc(vertexIndex, scene, hit_points[i], hit_normals[i], bounces + 1, sample_count, order, sh_buffer);
        //    }
        //}
    };

    // Get the total number of coefficients for a function represented by
// all spherical harmonic basis of degree <= @order (it is a point of
// confusion that the order of an SH refers to its degree and not the order).
    constexpr int GetCoefficientCount(int order) {
        return (order + 1) * (order + 1);
    }

    virtual int preprocess(const Scene *scene) override
    {

        // Here only compute one mesh
        const auto mesh = scene->getMeshes()[0];
        // Projection environment
        auto cubePath = getFileResolver()->resolve(m_CubemapPath);
        auto lightPath = cubePath / "light.txt";
        auto transPath = cubePath / "transport.txt";
        std::ofstream lightFout(lightPath.str());
        std::ofstream fout(transPath.str());
        int width, height, channel;
        std::vector<std::unique_ptr<float[]>> images =
            ProjEnv::LoadCubemapImages(cubePath.str(), width, height, channel);
        auto envCoeffs = ProjEnv::PrecomputeCubemapSH<SHOrder>(images, width, height, channel);
        m_LightCoeffs.resize(3, SHCoeffLength);
        for (int i = 0; i < envCoeffs.size(); i++)
        {
            lightFout << (envCoeffs)[i].x() << " " << (envCoeffs)[i].y() << " " << (envCoeffs)[i].z() << std::endl;
            m_LightCoeffs.col(i) = (envCoeffs)[i];
        }
        std::cout << "Computed light sh coeffs from: " << cubePath.str() << " to: " << lightPath.str() << std::endl;
        // Projection transport
        m_TransportSHCoeffs.resize(SHCoeffLength, mesh->getVertexCount());
        fout << mesh->getVertexCount() << std::endl;
        const double PI = 3.1415926;

        
        for (int i = 0; i < mesh->getVertexCount(); i++)
        {
            const Point3f &v = mesh->getVertexPositions().col(i);
            const Normal3f &n = mesh->getVertexNormals().col(i);

            auto shFunc = [&](double phi, double theta) -> double {
                Eigen::Array3d d = sh::ToVector(phi, theta);
                const auto wi = Vector3f(d.x(), d.y(), d.z());
                if (m_Type == Type::Unshadowed)
                {
                    // TODO: here you need to calculate unshadowed transport term of a given direction
                    // TODO: 此处你需要计算给定方向下的unshadowed传输项球谐函数值
                    //s = std::max(0.0, 5 * cos(theta) - 4) + std::max(0.0, -4 * sin(theta - PI) * cos(phi - 2.5) - 3);

                    double s = 0;
                    float m = 0.0f;
                    s = std::max(0.0, (double)wi.dot(n));
                    //auto value = H * preGenVector[i].SH_value[j];
                    //s = 0.2;
                    //if (s > 0.0001) {
                    //    auto vs = 0;
                    //}
                    return s;
                }
                else
                {
                    // TODO: here you need to calculate shadowed transport term of a given direction
                    // TODO: 此处你需要计算给定方向下的shadowed传输项球谐函数值

                    double s = 0;
                    float m = 0.0f;
                    m = n.dot(wi);
                    double p = (double)m;
                    s = std::max(0.0, p);
                    Ray3f ray = Ray3f(v, wi);
                    bool cross = scene->rayIntersect(ray);
                    if (cross) {
                        s = 0;
                    }

                    //if (s > 0.0001) 
                    //{
                    //    auto vs = 0;
                    //}

                    return s;
                }
            };
            auto shCoeff = sh::ProjectFunction(SHOrder, shFunc, m_SampleCount);
            for (int j = 0; j < shCoeff->size(); j++)
            {
                m_TransportSHCoeffs.col(i).coeffRef(j) = (*shCoeff)[j];
            }
        }
        
        if (m_Type == Type::Interreflection)
        {
            // TODO: leave for bonus
            std::vector<Eigen::MatrixXf> sh_buffer;
            sh_buffer.reserve(m_Bounce + 1);
            sh_buffer.push_back(m_TransportSHCoeffs);
            for (int i = 1; i < m_Bounce + 1; i++)
            {
                Eigen::MatrixXf newMatrixXf;
                newMatrixXf.resize(SHCoeffLength, mesh->getVertexCount());
                sh_buffer.push_back(newMatrixXf);
            }

            std::vector<ProjEnv::InterRef> hits;
            for (int i = 0; i < mesh->getVertexCount(); i++)
            {                
                const Point3f& v = mesh->getVertexPositions().col(i);
                const Normal3f& n = mesh->getVertexNormals().col(i);

                interreflectionFunc(i,scene,v,n, 1, m_SampleCount, SHOrder, sh_buffer, hits);
            }
            for (int  b = 2; b <= m_Bounce; b++)
            {
                std::vector<ProjEnv::InterRef> newhits;
                for (int ref = 0; ref < hits.size(); ref++) {
                    interreflectionFunc(hits[ref].vertexIndex, scene, hits[ref].p, hits[ref].n, b, m_SampleCount, SHOrder, sh_buffer, newhits);
                }
                hits = newhits;
            }

            
            for (int i = 0; i < mesh->getVertexCount(); i++)
            {
                for (int b = 1; b < m_Bounce + 1; b++)
                {
                    for (int j = 0; j < sh_buffer[b].col(i).size(); j++)
                    {
                        m_TransportSHCoeffs.col(i).coeffRef(j) += sh_buffer[b].col(i).coeffRef(j);
                    }                     
                }
            }
            
        }

        // Save in face format
        for (int f = 0; f < mesh->getTriangleCount(); f++)
        {
            const MatrixXu &F = mesh->getIndices();
            uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);
            for (int j = 0; j < SHCoeffLength; j++)
            {
                fout << m_TransportSHCoeffs.col(idx0).coeff(j) << " ";
            }
            fout << std::endl;
            for (int j = 0; j < SHCoeffLength; j++)
            {
                fout << m_TransportSHCoeffs.col(idx1).coeff(j) << " ";
            }
            fout << std::endl;
            for (int j = 0; j < SHCoeffLength; j++)
            {
                fout << m_TransportSHCoeffs.col(idx2).coeff(j) << " ";
            }
            fout << std::endl;
        }
        std::cout << "Computed SH coeffs"
                  << " to: " << transPath.str() << std::endl;
        return 0;
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
    {
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        const Eigen::Matrix<Vector3f::Scalar, SHCoeffLength, 1> sh0 = m_TransportSHCoeffs.col(its.tri_index.x()),
                                                                sh1 = m_TransportSHCoeffs.col(its.tri_index.y()),
                                                                sh2 = m_TransportSHCoeffs.col(its.tri_index.z());
        const Eigen::Matrix<Vector3f::Scalar, SHCoeffLength, 1> rL = m_LightCoeffs.row(0), gL = m_LightCoeffs.row(1), bL = m_LightCoeffs.row(2);

        Color3f c0 = Color3f(rL.dot(sh0), gL.dot(sh0), bL.dot(sh0)),
                c1 = Color3f(rL.dot(sh1), gL.dot(sh1), bL.dot(sh1)),
                c2 = Color3f(rL.dot(sh2), gL.dot(sh2), bL.dot(sh2));

        const Vector3f &bary = its.bary;
        Color3f c = bary.x() * c0 + bary.y() * c1 + bary.z() * c2;
        // TODO: you need to delete the following four line codes after finishing your calculation to SH,
        //       we use it to visualize the normals of model for debug.
        // TODO: 在完成了球谐系数计算后，你需要删除下列四行，这四行代码的作用是用来可视化模型法线
        //if (c.isZero()) {
        //    auto n_ = its.shFrame.n.cwiseAbs();
        //    return Color3f(n_.x(), n_.y(), n_.z());
        //}
        return c;
    }

    std::string toString() const
    {
        return "PRTIntegrator[]";
    }

private:
    Type m_Type;
    int m_Bounce = 1;
    int m_SampleCount = 100;
    std::string m_CubemapPath;
    Eigen::MatrixXf m_TransportSHCoeffs;
    Eigen::MatrixXf m_LightCoeffs;
};

NORI_REGISTER_CLASS(PRTIntegrator, "prt");
NORI_NAMESPACE_END