#include <glm.hpp>
#include <ProgressBar.hpp>

#include "cyHairFile.h"
#include <array>
#include <exception>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

/****************************************************************************************************/
using Vec3f  = glm::tvec3<float>;
using Vec2ui = glm::tvec2<unsigned>;

unsigned              nStrands;
unsigned              nTotalVertices;
std::vector<Vec3f>    vertices;
std::vector<float>    thicknesses;
std::vector<Vec2ui>   strandRanges;      /* range of vertices in each strand */
std::vector<unsigned> strandNumSegments; /* number of segments in each strand */

Vec3f bMin, bMax, diff, center;

constexpr float thicknessScale  = 1.0f;
constexpr float thicknessJitter = 0.0f;
constexpr float yScale          = 1.0f;

constexpr float scale       = 0.17f;
constexpr Vec3f translation = Vec3f(0, 7.5, 0);

/****************************************************************************************************/
void loadCYHairModel(const std::string hairFile) {
    cyHairFile cyHair;
    if(cyHair.LoadFromFile(hairFile.c_str()) < 0) {
        throw  std::exception("Cannot load file!");
    }
    nStrands       = cyHair.GetHeader().hair_count;
    nTotalVertices = cyHair.GetHeader().point_count;

    vertices.resize(nTotalVertices);
    thicknesses.resize(nTotalVertices);
    ////////////////////////////////////////////////////////////////////////////////
    // populate data
    if(cyHair.GetPointsArray() == nullptr) {
        throw  std::exception("Cannot load point array!");
    }
    std::memcpy((void*)vertices.data(), (void*)cyHair.GetPointsArray(), nTotalVertices * sizeof(Vec3f));

    // swap y-z
    bMin = Vec3f(1e10);
    bMax = Vec3f(-1e10);
    for(auto& vertex : vertices) {
        std::swap(vertex.y, vertex.z);
        //        std::swap(vertex.x, vertex.z);
        vertex *= scale;
        vertex += translation;

        for(unsigned i = 0; i < 3; ++i) {
            if(bMin[i] > vertex[i]) {
                bMin[i] = vertex[i];
            }
            if(bMax[i] < vertex[i]) {
                bMax[i] = vertex[i];
            }
        }
    }
    diff   = bMax - bMin;
    center = (bMin + bMax) * 0.5f;

    {
        for(auto& vertex : vertices) {
            vertex[1] -= center.y;
            vertex[1] *= yScale;
            vertex[1] += center.y;
        }
    }

    bool bHasSegmentArray   = cyHair.GetSegmentsArray() != nullptr;
    bool bHasThicknessArray = cyHair.GetThicknessArray() != nullptr;
    strandNumSegments.resize(nStrands);
    strandRanges.resize(nStrands);

    auto rand11 = []() {
                      const auto tmp = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                      return tmp * 2 - 1;
                  };
    srand(time(nullptr));

    unsigned startVertex = 0;
    unsigned endVertex   = 0;
    for(unsigned idx = 0; idx < nStrands; ++idx) {
        strandNumSegments[idx] = bHasSegmentArray ?
                                 static_cast<unsigned>(cyHair.GetSegmentsArray()[idx]) :
                                 cyHair.GetHeader().d_segments;
        endVertex         = startVertex + strandNumSegments[idx] + 1;
        strandRanges[idx] = Vec2ui{ startVertex, endVertex };

        for(unsigned j = startVertex; j < endVertex; ++j) {
            thicknesses[j] = bHasThicknessArray ?
                             cyHair.GetThicknessArray()[j] :
                             cyHair.GetHeader().d_thickness;
            thicknesses[j] *= (scale * thicknessScale);
            if(thicknessJitter > 0) {
                thicknesses[j] *= (1 + thicknessJitter * rand11());
            }
        }
        startVertex = endVertex;
    }

    std::cout << "Loaded file " << hairFile << std::endl;
}

void savePBRT(const std::string fileName) {
    std::ofstream file(fileName);
    if(!file.is_open()) {
        throw std::exception("Cannot open file for writing!");
    }

    std::cout << "Converting file " << fileName << std::endl;
    using VPoints = std::array<Vec3f, 4>;
    const char* curveLine = R"(Shape "curve" "string type" [ "cylinder" ] "point P" [ %f %f %f %f %f %f %f %f %f %f %f %f ] "float width0"  [ %f ] "float width1" [ %f ])";

    auto writeFile = [&](const VPoints& points, float t0, float t1) {
                         char buff[256];
                         sprintf_s(buff,
                                   curveLine,
                                   points[0].x,
                                   points[0].y,
                                   points[0].z,
                                   //
                                   points[1].x,
                                   points[1].y,
                                   points[1].z,
                                   //
                                   points[2].x,
                                   points[2].y,
                                   points[2].z,
                                   //
                                   points[3].x,
                                   points[3].y,
                                   points[3].z,
                                   //
                                   t0, t1
                                   );
                         file << buff;
                         file << "\n";
                     };

    auto convert =
        [](const VPoints& controlPoints) {
            static constexpr float alpha = 0.5f;

            const auto P0 = controlPoints[0];
            const auto P1 = controlPoints[0 + 1];
            const auto P2 = controlPoints[0 + 2];
            const auto P3 = controlPoints[0 + 3];

            const float d1 = (P0 - P1).length();
            const float d2 = (P1 - P2).length();
            const float d3 = (P2 - P3).length();

            const float d1_alpha  = std::pow(d1, alpha);
            const float d2_alpha  = std::pow(d2, alpha);
            const float d3_alpha  = std::pow(d3, alpha);
            const float d1_2alpha = std::pow(d1, 2.0f * alpha);
            const float d2_2alpha = std::pow(d2, 2.0f * alpha);
            const float d3_2alpha = std::pow(d3, 2.0f * alpha);

            /* Compute control points of the cubic Bezier curve */
            VPoints Bi;
            Bi[0] = P1;
            Bi[3] = P2;
            Bi[1] = (d1_2alpha * P2 - d2_2alpha * P0 + (2.0f * d1_2alpha + 3.0f * d1_alpha * d2_alpha + d2_2alpha) * P1) /
                    (3.0f * d1_alpha * (d1_alpha + d2_alpha));
            Bi[2] = (d3_2alpha * P1 - d2_2alpha * P3 + (2.0f * d3_2alpha + 3.0f * d3_alpha * d2_alpha + d2_2alpha) * P2) /
                    (3.0f * d3_alpha * (d3_alpha + d2_alpha));
            return Bi;
        };

    progresscpp::ProgressBar progressBar(nStrands, 70, '#', '-');
    unsigned                 count = 0;
    for(const auto& range: strandRanges) {
        ++progressBar;
        if(count % (nStrands / 10u) == 0) {
            progressBar.display();
        }
        if(range[1] - range[0] < 4) {
            continue;
        }

        // First segment
        {
            const VPoints points = { vertices[range[0]],
                                     vertices[range[0]],
                                     vertices[range[0] + 1u],
                                     vertices[range[0] + 2u] };
            const auto    cubicControlPoints = convert(points);
            writeFile(cubicControlPoints, thicknesses[range[0]], thicknesses[range[0] + 2u]);
        }

        for(unsigned i = range[0] + 1u; i < range[1] - 3u; ++i) {
            const VPoints points = { vertices[i],
                                     vertices[i + 1u],
                                     vertices[i + 2u],
                                     vertices[i + 3u] };
            const auto    cubicControlPoints = convert(points);
            writeFile(cubicControlPoints, thicknesses[i], thicknesses[i + 3u]);
        }

        // Last segment
        {
            const VPoints points = { vertices[range[1] - 3u],
                                     vertices[range[1] - 2u],
                                     vertices[range[1] - 1u],
                                     vertices[range[1] - 1u] };
            const auto    cubicControlPoints = convert(points);
            writeFile(cubicControlPoints, thicknesses[range[1] - 3u], thicknesses[range[1] - 1u]);
        }
    }

    file.close();

    std::cout << "Process " << nStrands << std::endl;
    std::cout << "Vertices: " << nTotalVertices << std::endl;
    std::cout << "Min: " << bMin.x << ", " << bMin.y << ", " << bMin.z << std::endl;
    std::cout << "Max: " << bMax.x << ", " << bMax.y << ", " << bMax.z << std::endl;
    std::cout << "Diff: " << diff.x << ", " << diff.y << ", " << diff.z << std::endl;
    std::cout << "Center: " << center.x << ", " << center.y << ", " << center.z << std::endl;
}

int main() {
#if 0
    const std::string fileName = "wWavy";
    /*
       Process 50000
       Vertices: 2450000
       Min: -60.8702, -56.9884, -46.3591
       Max: 31.4414, 66.5383, 45.5207
       Diff: 92.3116, 123.527, 91.8797
       Center: -14.7144, 4.77494, -0.419195

       constexpr float thicknessScale  = 1.0f;
       constexpr float thicknessJitter = 0.0f;
       constexpr float yScale          = 1.0f;

       constexpr float scale       = 0.13f;
       constexpr Vec3f translation = Vec3f(0, 9.5, 0);

       Min: -6.02668, 2.09151, -7.91313
       Max: 5.91769, 18.15, 4.08738
       Diff: 11.9444, 16.0585, 12.0005
       Center: -0.0544953, 10.1207, -1.91287
     */
#endif

#if 0
    const std::string fileName = "wWavyThin";
    /*
       Process 10000
       Vertices: 872756
       Min: -72.3699, -37.8153, -44.5626
       Max: 22.8397, 67.7453, 43.7826
       Diff: 95.2096, 105.561, 88.3451
       Center: -24.7651, 14.965, -0.390005

       constexpr float thicknessScale  = 1.0f;
       constexpr float thicknessJitter = 0.0f;
       constexpr float yScale          = 1.0f;

       constexpr float scale       = 0.13f;
       constexpr Vec3f translation = Vec3f(0, 9.5, 0);

       Min: -5.79313, 4.58401, -9.40809
       Max: 5.69173, 18.3069, 2.96916
       Diff: 11.4849, 13.7229, 12.3772
       Center: -0.0507007, 11.4455, -3.21946
     */
#endif

#if 1
    const std::string fileName = "dark";
    /*
       Process 15000
       Vertices: 1031268
       Min: -35.3979, -23.4088, -43.6086
       Max: 33.6512, 63.4912, 19.1387
       Diff: 69.0491, 86.9, 62.7474
       Center: -0.873335, 20.0412, -12.235

       No swap x-z

       constexpr float thicknessScale  = 1.0f;
       constexpr float thicknessJitter = 0.0f;
       constexpr float yScale          = 1.0f;

       constexpr float scale       = 0.17f;
       constexpr Vec3f translation = Vec3f(0, 7.5, 0);

       Min: -6.01764, 3.52051, -7.41347
       Max: 5.7207, 18.2935, 3.25358
       Diff: 11.7383, 14.773, 10.6671
       Center: -0.148467, 10.907, -2.07994
     */
#endif

#if 0
    const std::string fileName = "long-straight-hair";
    /*
       Process 10000
       Vertices: 160000
       Min: -32.4956, -22.7086, -33.9009
       Max: 30.8987, 63.678, 24.074
       Diff: 63.3943, 86.3865, 57.9749
       Center: -0.798452, 20.4847, -4.91345

       constexpr float thicknessScale  = 1.0f;
       constexpr float thicknessJitter = 0.0f;
       constexpr float yScale          = 1.0f;

       constexpr float scale       = 0.18f;
       constexpr Vec3f translation = Vec3f(0, 7.0, 0);

       Min: -5.84921, 2.91246, -6.10216
       Max: 5.56177, 18.462, 4.33332
       Diff: 11.411, 15.5496, 10.4355
       Center: -0.143722, 10.6872, -0.884421
     */
#endif

#if 0
    const std::string fileName = "natural";
    /*

       constexpr float thicknessScale  = 1.0f;
       constexpr float thicknessJitter = 0.0f;
       constexpr float yScale          = 1.0f;

       constexpr float scale       = 0.15f;
       constexpr Vec3f translation = Vec3f(0, 8, 0);

       No x-z swap

       Process 10000
       Vertices: 1519823
       Min: -6.89381, 4.14137, -9.13453
       Max: 6.39393, 17.8792, 3.25541
       Diff: 13.2877, 13.7378, 12.3899
       Center: -0.249938, 11.0103, -2.93956
     */
#endif

#if 0
    const std::string fileName = "wCurly";
    /*
       constexpr float scale       = 0.12f;
       constexpr Vec3f translation = Vec3f(0, 8, 0);


       constexpr float scale       = 0.1f;
       constexpr Vec3f translation = Vec3f(0, 10, 0);

       Process 50000
       Vertices: 3441580
       Min: -5.89147, 3.66136, -9.28943
       Max: 6.26518, 17.9824, 3.17144
       Diff: 12.1566, 14.3211, 12.4609
       Center: 0.186856, 10.8219, -3.059
     */
#endif

#if 0
    const std::string fileName = "wStraight";
    /*
       constexpr float thicknessScale  = 1.0f;
       constexpr float thicknessJitter = 0.0f;
       constexpr float yScale          = 1.0f;

       constexpr float scale       = 0.13f;
       constexpr Vec3f translation = Vec3f(0, 9.5, 0);

       Process 50000
       Vertices: 3441580
       Min: -4.25692, 4.28597, -5.85832
       Max: 4.33566, 15.9567, 2.85001
       Diff: 8.59258, 11.6707, 8.70833
       Center: 0.0393696, 10.1213, -1.50416
     */
#endif

    /****************************************************************************************************/
#if 0
    const std::string fileName = "wStraight";
    /*
       constexpr float scale          = 0.14f;
       constexpr float thicknessScale = 0.35f;
       constexpr float yScale         = 0.7f;
       constexpr Vec3f translation    = Vec3f(0, 10, 0);

       Process 50000
       Vertices: 1250000
       Min: -5.95969, 2.70035, -8.20165
       Max: 6.06992, 19.0393, 3.99001
       Diff: 12.0296, 16.339, 12.1917
       Center: 0.0551174, 10.8698, -2.10582
     */
#endif

#if 0
    const std::string fileName = "wCurly";
    /*
       constexpr float scale          = 0.12f;
       constexpr float thicknessScale = 0.35f;
       constexpr float yScale         = 0.7f;
       constexpr Vec3f translation    = Vec3f(0, 10, 0);

       Process 50000
       Vertices: 3441580
       Min: -5.89147, 3.66136, -9.28943
       Max: 6.26518, 17.9824, 3.17144
       Diff: 12.1566, 14.3211, 12.4609
       Center: 0.186856, 10.8219, -3.059
     */
#endif

#if 0
    const std::string fileName = "long-straight-hair";
    /*
       Process 10000
       Vertices: 160000
       Min: -32.4956, -22.7086, -33.9009
       Max: 30.8987, 63.678, 24.074
       Diff: 63.3943, 86.3865, 57.9749
       Center: -0.798452, 20.4847, -4.91345

       constexpr float scale          = 0.19f;
       constexpr float thicknessScale = 0.5f;
       constexpr float yScale         = 0.7f;
       constexpr Vec3f translation    = Vec3f(0, 7, 0);

       Min: -4.54938, 4.3208, -4.74612
       Max: 4.32582, 16.4149, 3.37036
       Diff: 8.8752, 12.0941, 8.11648
       Center: -0.111783, 10.3679, -0.687883
     */
#endif

#if 0
    const std::string fileName = "natural";
    /*
       constexpr float scale          = 0.153f;
       constexpr float thicknessScale = 0.35f;
       constexpr float yScale         = 0.7f;
       constexpr Vec3f translation    = Vec3f(0, 8, 0);

       No x-z swap

       Process 10000
       Vertices: 1519823
       Min: -7.03169, 4.06419, -9.31722
       Max: 6.52181, 18.0767, 3.32052
       Diff: 13.5535, 14.0126, 12.6377
       Center: -0.254937, 11.0705, -2.99835
     */
#endif

    loadCYHairModel(fileName + ".hair");
    savePBRT(fileName + ".pbrt");
    return 0;
}
