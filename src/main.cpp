#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

namespace Physics
{
    constexpr double G = 6.67430e-11;
    constexpr double SUN_MASS = 1.98847e30;
    constexpr double EARTH_MASS = 5.9722e24;
    constexpr double AU_M = 149'597'870'700.0;
    constexpr double AU_KM = AU_M / 1000.0;
    constexpr double DAY_S = 86'400.0;
    constexpr double PI = 3.1415926535897932384626433832795;

    constexpr double ADITYA_L1_REFERENCE_KM = 1'500'000.0;
    constexpr double ADITYA_HALO_PERIOD_DAYS = 177.86;
    constexpr double MISSION_LIFETIME_YEARS = 5.0;

    struct LPoint
    {
        std::string name;
        double x = 0.0;
        double y = 0.0;
        double earthKm = 0.0;
        double sunKm = 0.0;
    };

    struct State
    {
        glm::dvec3 r{};
        glm::dvec3 v{};
    };
}

static double Mu()
{
    return Physics::EARTH_MASS /
           (Physics::SUN_MASS + Physics::EARTH_MASS);
}

static double LagrangeEquation(double x)
{
    const double mu = Mu();
    const double sunX = -mu;
    const double earthX = 1.0 - mu;

    const double r1 = std::abs(x - sunX);
    const double r2 = std::abs(x - earthX);

    if (r1 < 1.0e-14 || r2 < 1.0e-14)
        return std::numeric_limits<double>::quiet_NaN();

    return x
        - (1.0 - mu) * (x - sunX) / (r1 * r1 * r1)
        - mu * (x - earthX) / (r2 * r2 * r2);
}

static double SolveRootBracketed(double lo, double hi)
{
    double flo = LagrangeEquation(lo);
    double fhi = LagrangeEquation(hi);

    if (!std::isfinite(flo) || !std::isfinite(fhi))
        return 0.0;

    if (flo == 0.0)
        return lo;
    if (fhi == 0.0)
        return hi;

    if (flo * fhi > 0.0)
    {
        // Deterministic scan inside the interval. This avoids startup
        // failures when a very narrow Sun-Earth root is encountered.
        constexpr int samples = 16384;
        double px = lo;
        double pf = flo;

        for (int i = 1; i <= samples; ++i)
        {
            const double x = lo + (hi - lo) * static_cast<double>(i) / samples;
            const double fx = LagrangeEquation(x);
            if (std::isfinite(fx) && pf * fx <= 0.0)
            {
                lo = px;
                hi = x;
                flo = pf;
                fhi = fx;
                break;
            }
            px = x;
            pf = fx;
        }
    }

    if (flo * fhi > 0.0)
    {
        // Newton fallback from the midpoint; bounded to the requested interval.
        double x = 0.5 * (lo + hi);
        for (int i = 0; i < 200; ++i)
        {
            const double f = LagrangeEquation(x);
            const double h = 1.0e-7;
            const double fp = (LagrangeEquation(x + h) - LagrangeEquation(x - h)) / (2.0 * h);
            if (!std::isfinite(f) || !std::isfinite(fp) || std::abs(fp) < 1.0e-15)
                break;
            const double next = x - f / fp;
            if (!std::isfinite(next) || next < lo || next > hi)
                break;
            x = next;
            if (std::abs(f) < 1.0e-13)
                return x;
        }

        // Known closed-form geometry is not needed for the collinear points;
        // this deterministic seed is only a final guard against a fatal startup.
        return 0.5 * (lo + hi);
    }

    for (int i = 0; i < 300; ++i)
    {
        const double mid = 0.5 * (lo + hi);
        const double fm = LagrangeEquation(mid);

        if (!std::isfinite(fm))
            break;

        if (std::abs(fm) < 1.0e-14)
            return mid;

        if (flo * fm <= 0.0)
        {
            hi = mid;
            fhi = fm;
        }
        else
        {
            lo = mid;
            flo = fm;
        }
    }

    return 0.5 * (lo + hi);
}

static std::array<Physics::LPoint, 5> ComputeLagrangePoints()
{
    const double mu = Mu();
    const double earthX = 1.0 - mu;

    std::array<Physics::LPoint, 5> p{};

    // Direct non-singular brackets for the Sun-Earth CR3BP.
    p[0] = {"L1", SolveRootBracketed(earthX - 0.10, earthX - 1.0e-8), 0.0, 0.0, 0.0};
    p[1] = {"L2", SolveRootBracketed(earthX + 1.0e-8, earthX + 0.10), 0.0, 0.0, 0.0};
    p[2] = {"L3", SolveRootBracketed(-1.50, -1.0), 0.0, 0.0, 0.0};
    p[3] = {"L4", 0.5 - mu, std::sqrt(3.0) / 2.0, 0.0, 0.0};
    p[4] = {"L5", 0.5 - mu, -std::sqrt(3.0) / 2.0, 0.0, 0.0};

    for (auto& q : p)
    {
        q.earthKm = std::hypot(q.x - earthX, q.y) * Physics::AU_KM;
        q.sunKm = std::hypot(q.x + mu, q.y) * Physics::AU_KM;
    }

    return p;
}

static double MeanMotion()
{
    return std::sqrt(
        Physics::G * (Physics::SUN_MASS + Physics::EARTH_MASS) /
        std::pow(Physics::AU_M, 3.0)
    );
}

static double TimeScaleDays()
{
    return (1.0 / MeanMotion()) / Physics::DAY_S;
}

static double VelocityScaleKmS()
{
    return Physics::AU_KM * MeanMotion();
}

// Exact time-dependent inertial visualization for the circular Sun-Earth
// CR3BP. The L1-L5 coordinates are fixed in the rotating frame; in an
// inertial frame the Sun-Earth system rotates with mean motion n.
static double OrbitalAngleRad(double simDays)
{
    return MeanMotion() * simDays * Physics::DAY_S;
}

// Forward declaration: display mapping is defined below.
static glm::vec3 DisplayPoint(const glm::dvec3& n);

static glm::vec3 RotateDisplay(const glm::vec3& p, double angle)
{
    const float c = static_cast<float>(std::cos(angle));
    const float s = static_cast<float>(std::sin(angle));
    return {
        c * p.x - s * p.y,
        s * p.x + c * p.y,
        p.z
    };
}

static glm::vec3 InertialDisplayPoint(const glm::dvec3& n, double simDays)
{
    return RotateDisplay(DisplayPoint(n), OrbitalAngleRad(simDays));
}

static glm::dvec3 RotateNormalized(const glm::dvec3& n, double angle)
{
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return {
        c * n.x - s * n.y,
        s * n.x + c * n.y,
        n.z
    };
}

static glm::vec3 CurrentSunDisplay(double simDays)
{
    return InertialDisplayPoint({-Mu(), 0.0, 0.0}, simDays);
}

static glm::vec3 CurrentEarthDisplay(double simDays)
{
    return InertialDisplayPoint({1.0 - Mu(), 0.0, 0.0}, simDays);
}

static double EarthDistanceKm(const glm::dvec3& n)
{
    const double earthX = 1.0 - Mu();
    return std::hypot(n.x - earthX, n.y) * Physics::AU_KM;
}

static double SunDistanceKm(const glm::dvec3& n)
{
    return std::hypot(n.x + Mu(), n.y) * Physics::AU_KM;
}

// Display mapping: 1 AU maps to 3 display units.
// Physics values remain in normalized CR3BP coordinates.
static glm::vec3 DisplayPoint(const glm::dvec3& n)
{
    constexpr float scale = 3.0f;
    return {
        static_cast<float>((n.x - 0.5) * scale),
        static_cast<float>(n.y * scale),
        static_cast<float>(n.z * scale)
    };
}

static glm::dvec3 NormalizedPoint(const Physics::LPoint& p)
{
    return {p.x, p.y, 0.0};
}

static glm::dvec3 CR3BPAcceleration(const Physics::State& s)
{
    const double mu = Mu();
    const double x = s.r.x;
    const double y = s.r.y;
    const double z = s.r.z;
    const double vx = s.v.x;
    const double vy = s.v.y;

    const double r1sq = (x + mu) * (x + mu) + y * y + z * z;
    const double r2sq = (x - (1.0 - mu)) * (x - (1.0 - mu)) + y * y + z * z;
    const double r1 = std::sqrt(r1sq);
    const double r2 = std::sqrt(r2sq);

    const double r13 = r1sq * r1;
    const double r23 = r2sq * r2;

    return {
        x - (1.0 - mu) * (x + mu) / r13 - mu * (x - (1.0 - mu)) / r23 + 2.0 * vy,
        y - (1.0 - mu) * y / r13 - mu * y / r23 - 2.0 * vx,
        -(1.0 - mu) * z / r13 - mu * z / r23
    };
}

static Physics::State RK4Step(const Physics::State& s, double dt)
{
    const auto f = [](const Physics::State& a)
    {
        return std::pair<glm::dvec3, glm::dvec3>{a.v, CR3BPAcceleration(a)};
    };

    const auto k1 = f(s);
    Physics::State s2{s.r + 0.5 * dt * k1.first, s.v + 0.5 * dt * k1.second};
    const auto k2 = f(s2);
    Physics::State s3{s.r + 0.5 * dt * k2.first, s.v + 0.5 * dt * k2.second};
    const auto k3 = f(s3);
    Physics::State s4{s.r + dt * k3.first, s.v + dt * k3.second};
    const auto k4 = f(s4);

    return {
        s.r + (dt / 6.0) * (k1.first + 2.0 * k2.first + 2.0 * k3.first + k4.first),
        s.v + (dt / 6.0) * (k1.second + 2.0 * k2.second + 2.0 * k3.second + k4.second)
    };
}

static double Jacobi(const Physics::State& s)
{
    const double mu = Mu();
    const double r1 = std::sqrt((s.r.x + mu) * (s.r.x + mu) + s.r.y * s.r.y + s.r.z * s.r.z);
    const double r2 = std::sqrt((s.r.x - (1.0 - mu)) * (s.r.x - (1.0 - mu)) + s.r.y * s.r.y + s.r.z * s.r.z);
    const double omega = 0.5 * (s.r.x * s.r.x + s.r.y * s.r.y) + (1.0 - mu) / r1 + mu / r2;
    return 2.0 * omega - glm::dot(s.v, s.v);
}

struct Mesh
{
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint ebo = 0;
    GLsizei count = 0;
};

struct LineBuffer
{
    GLuint vao = 0;
    GLuint vbo = 0;
    std::size_t capacity = 0;
};

static GLuint CompileShader(GLenum type, const char* source)
{
    const GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);

    GLint ok = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok)
    {
        char log[4096]{};
        glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
        std::cerr << "Shader error: " << log << '\n';
        glDeleteShader(shader);
        return 0;
    }
    return shader;
}

static GLuint CreateProgram(const char* vs, const char* fs)
{
    const GLuint v = CompileShader(GL_VERTEX_SHADER, vs);
    const GLuint f = CompileShader(GL_FRAGMENT_SHADER, fs);
    if (!v || !f)
        return 0;

    const GLuint p = glCreateProgram();
    glAttachShader(p, v);
    glAttachShader(p, f);
    glLinkProgram(p);

    GLint ok = 0;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok)
    {
        char log[4096]{};
        glGetProgramInfoLog(p, sizeof(log), nullptr, log);
        std::cerr << "Program link error: " << log << '\n';
        glDeleteProgram(p);
        glDeleteShader(v);
        glDeleteShader(f);
        return 0;
    }

    glDeleteShader(v);
    glDeleteShader(f);
    return p;
}

static Mesh MakeSphere(int stacks, int slices)
{
    struct Vertex { glm::vec3 p; glm::vec3 n; };
    std::vector<Vertex> vertices;
    std::vector<uint32_t> indices;
    vertices.reserve(static_cast<std::size_t>((stacks + 1) * (slices + 1)));

    for (int y = 0; y <= stacks; ++y)
    {
        const float v = static_cast<float>(y) / static_cast<float>(stacks);
        const float phi = v * static_cast<float>(Physics::PI);
        for (int x = 0; x <= slices; ++x)
        {
            const float u = static_cast<float>(x) / static_cast<float>(slices);
            const float theta = u * 2.0f * static_cast<float>(Physics::PI);
            const glm::vec3 n{
                std::sin(phi) * std::cos(theta),
                std::cos(phi),
                std::sin(phi) * std::sin(theta)
            };
            vertices.push_back({n, n});
        }
    }

    for (int y = 0; y < stacks; ++y)
        for (int x = 0; x < slices; ++x)
        {
            const uint32_t a = static_cast<uint32_t>(y * (slices + 1) + x);
            const uint32_t b = a + static_cast<uint32_t>(slices + 1);
            indices.insert(indices.end(), {a, b, a + 1, b, b + 1, a + 1});
        }

    Mesh mesh;
    glGenVertexArrays(1, &mesh.vao);
    glGenBuffers(1, &mesh.vbo);
    glGenBuffers(1, &mesh.ebo);

    glBindVertexArray(mesh.vao);
    glBindBuffer(GL_ARRAY_BUFFER, mesh.vbo);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(vertices.size() * sizeof(Vertex)), vertices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(indices.size() * sizeof(uint32_t)), indices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, n)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    mesh.count = static_cast<GLsizei>(indices.size());
    return mesh;
}

static Mesh MakeCube()
{
    const float s = 0.5f;
    const float data[] = {
        -s,-s,-s, 0,0,-1,  s,-s,-s, 0,0,-1,  s,s,-s,0,0,-1,  -s,s,-s,0,0,-1,
        -s,-s,s, 0,0,1,  s,-s,s,0,0,1,  s,s,s,0,0,1,  -s,s,s,0,0,1,
        -s,-s,-s,-1,0,0, -s,s,-s,-1,0,0, -s,s,s,-1,0,0, -s,-s,s,-1,0,0,
        s,-s,-s,1,0,0, s,s,-s,1,0,0, s,s,s,1,0,0, s,-s,s,1,0,0,
        -s,-s,-s,0,-1,0, -s,-s,s,0,-1,0, s,-s,s,0,-1,0, s,-s,-s,0,-1,0,
        -s,s,-s,0,1,0, -s,s,s,0,1,0, s,s,s,0,1,0, s,s,-s,0,1,0
    };
    const uint32_t indices[] = {
        0,1,2, 2,3,0, 4,6,5, 6,4,7, 8,9,10, 10,11,8,
        12,14,13, 14,12,15, 16,17,18, 18,19,16, 20,22,21, 22,20,23
    };

    Mesh mesh;
    glGenVertexArrays(1, &mesh.vao);
    glGenBuffers(1, &mesh.vbo);
    glGenBuffers(1, &mesh.ebo);
    glBindVertexArray(mesh.vao);
    glBindBuffer(GL_ARRAY_BUFFER, mesh.vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(data), data, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    mesh.count = 36;
    return mesh;
}

static void DestroyMesh(Mesh& mesh)
{
    if (mesh.ebo) glDeleteBuffers(1, &mesh.ebo);
    if (mesh.vbo) glDeleteBuffers(1, &mesh.vbo);
    if (mesh.vao) glDeleteVertexArrays(1, &mesh.vao);
    mesh = {};
}

static void DestroyLineBuffer(LineBuffer& b)
{
    if (b.vbo) glDeleteBuffers(1, &b.vbo);
    if (b.vao) glDeleteVertexArrays(1, &b.vao);
    b = {};
}

struct Renderer
{
    GLFWwindow* window = nullptr;
    GLuint meshProgram = 0;
    GLuint pointProgram = 0;
    GLuint lineProgram = 0;
    Mesh sphere{};
    Mesh spacecraftBody{};
    Mesh spacecraftPanel{};
    LineBuffer dynamicLines{};
    LineBuffer dynamicPoints{};
    std::array<Physics::LPoint, 5> points{};
    std::vector<glm::vec3> stars;
    std::vector<glm::vec3> trajectory;

    bool showPoints = true;
    bool showHalo = true;
    bool showTrajectory = false;
    bool showFlares = true;
    bool showRadiation = true;
    bool showSolarWind = true;
    bool showMagnetic = true;
    bool showSolarMagnetic = true;
    bool showInfo = true;
    bool showWorldLabels = true;
    bool followSpacecraft = false;
    bool paused = false;
    bool overview = true;
    bool experimental = false;

    int selectedPoint = 0;
    int positionMode = 0;
    glm::dvec3 experimentalPosition{0.0, 0.0, 0.0};

    float yaw = 0.0f;
    float pitch = 0.0f;
    float distance = 8.5f;
    bool rotatingCamera = false;
    double lastMouseX = 0.0;
    double lastMouseY = 0.0;

    double simDays = 0.0;
    double simDaysPerSecond = 0.15;
    double flareStrength = 0.95;
    double radiationStrength = 0.70;
    double windDensity = 1.0;
};

static const char* MeshVS = R"GLSL(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNormal;
uniform mat4 uModel;
uniform mat4 uView;
uniform mat4 uProj;
uniform float uTime;
struct VSOut { vec3 worldPos; vec3 normal; };
out vec3 vWorldPos;
out vec3 vNormal;
void main(){
    vec4 wp=uModel*vec4(aPos,1.0);
    vWorldPos=wp.xyz;
    vNormal=mat3(transpose(inverse(uModel)))*aNormal;
    gl_Position=uProj*uView*wp;
}
)GLSL";

static const char* MeshFS = R"GLSL(
#version 330 core
in vec3 vWorldPos;
in vec3 vNormal;
uniform vec3 uColor;
uniform vec3 uLightDir;
uniform float uEmission;
uniform float uRim;
uniform float uTime;
uniform int uBody;
out vec4 FragColor;
float noise3(vec3 p){
    return 0.5+0.5*sin(p.x*7.1+p.y*11.7+p.z*5.3+sin(p.x*3.7+p.z*4.1)*2.0);
}
void main(){
    vec3 N=normalize(vNormal);
    float diffuse=max(dot(N,normalize(uLightDir)),0.0);
    float rim=pow(1.0-max(dot(N,vec3(0,0,1)),0.0),2.0);
    vec3 c=uColor;
    if(uBody==1){
        float g=noise3(N*4.0)+0.45*noise3(N*12.0);
        c*=0.72+0.55*g;
        c+=vec3(0.20,0.055,0.0)*(0.35+0.65*g);
    } else if(uBody==2){
        float lat=asin(clamp(N.y,-1.0,1.0));
        float lon=atan(N.z,N.x);
        float land=sin(lon*4.5+sin(lat*4.0)*1.2)+sin(lon*9.0+lat*6.0);
        float cloud=noise3(N*18.0);
        vec3 ocean=vec3(0.02,0.12,0.48);
        vec3 green=vec3(0.06,0.32,0.07);
        vec3 ice=vec3(0.74,0.86,0.96);
        c=mix(ocean,green,smoothstep(0.55,1.15,land));
        c=mix(c,ice,smoothstep(0.82,1.0,abs(lat)));
        c+=vec3(0.82)*smoothstep(0.72,0.90,cloud)*0.18;
    }
    c=c*(0.16+0.84*diffuse)+c*uEmission+c*rim*uRim;
    FragColor=vec4(c,1.0);
}
)GLSL";

static const char* PointVS = R"GLSL(
#version 330 core
layout(location=0) in vec3 aPos;
uniform mat4 uView; uniform mat4 uProj; uniform float uSize;
void main(){gl_Position=uProj*uView*vec4(aPos,1.0);gl_PointSize=uSize;}
)GLSL";

static const char* PointFS = R"GLSL(
#version 330 core
uniform vec3 uColor; uniform float uAlpha;
out vec4 FragColor;
void main(){vec2 p=gl_PointCoord-vec2(0.5);float d=length(p);if(d>0.5) discard;float a=(1.0-smoothstep(0.02,0.5,d))*uAlpha;FragColor=vec4(uColor,a);}
)GLSL";

static const char* LineVS = R"GLSL(
#version 330 core
layout(location=0) in vec3 aPos;
uniform mat4 uView; uniform mat4 uProj;
void main(){gl_Position=uProj*uView*vec4(aPos,1.0);}
)GLSL";

static const char* LineFS = R"GLSL(
#version 330 core
uniform vec3 uColor; uniform float uAlpha;
out vec4 FragColor;
void main(){FragColor=vec4(uColor,uAlpha);}
)GLSL";

static glm::vec3 CameraPosition(const Renderer& r)
{
    return {
        r.distance * std::cos(r.pitch) * std::sin(r.yaw),
        r.distance * std::sin(r.pitch),
        r.distance * std::cos(r.pitch) * std::cos(r.yaw)
    };
}

static glm::vec3 MissionSpacecraftDisplay(const Renderer& r)
{
    if (r.experimental)
        return InertialDisplayPoint(r.experimentalPosition, r.simDays);

    const glm::dvec3 l1Rotating{r.points[0].x, r.points[0].y, 0.0};
    const double phase = 2.0 * Physics::PI * (r.simDays / Physics::ADITYA_HALO_PERIOD_DAYS);

    // Illustrative halo around the mathematically exact L1 location in
    // the rotating CR3BP frame. The whole configuration is then rotated
    // into the inertial display frame with the Earth-Sun mean motion.
    const double ax = 0.0016 * std::cos(phase);
    const double ay = 0.0042 * std::sin(phase);
    const double az = 0.0028 * std::sin(2.0 * phase);

    const glm::dvec3 q = l1Rotating + glm::dvec3(ax, ay, az);
    return InertialDisplayPoint(q, r.simDays);
}

static void DrawMesh(Renderer& r,const Mesh& mesh,const glm::vec3& position,const glm::vec3& scale,const glm::vec3& color,float emission,float rim,float rotation,int body,const glm::mat4& view,const glm::mat4& projection,double time)
{
    glUseProgram(r.meshProgram);
    glm::mat4 model(1.0f);
    model=glm::translate(model,position);
    model=glm::rotate(model,rotation,glm::vec3(0.0f,1.0f,0.0f));
    model=glm::scale(model,scale);

    glUniformMatrix4fv(glGetUniformLocation(r.meshProgram,"uModel"),1,GL_FALSE,glm::value_ptr(model));
    glUniformMatrix4fv(glGetUniformLocation(r.meshProgram,"uView"),1,GL_FALSE,glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(r.meshProgram,"uProj"),1,GL_FALSE,glm::value_ptr(projection));
    glUniform3fv(glGetUniformLocation(r.meshProgram,"uColor"),1,glm::value_ptr(color));
    glUniform3f(glGetUniformLocation(r.meshProgram,"uLightDir"),-0.7f,0.25f,0.6f);
    glUniform1f(glGetUniformLocation(r.meshProgram,"uEmission"),emission);
    glUniform1f(glGetUniformLocation(r.meshProgram,"uRim"),rim);
    glUniform1f(glGetUniformLocation(r.meshProgram,"uTime"),static_cast<float>(time));
    glUniform1i(glGetUniformLocation(r.meshProgram,"uBody"),body);

    glBindVertexArray(mesh.vao);
    glDrawElements(GL_TRIANGLES,mesh.count,GL_UNSIGNED_INT,nullptr);
    glBindVertexArray(0);
}

static void UploadBuffer(LineBuffer& b,const std::vector<glm::vec3>& data,GLenum usage)
{
    if(!b.vao){glGenVertexArrays(1,&b.vao);glGenBuffers(1,&b.vbo);}
    glBindVertexArray(b.vao);
    glBindBuffer(GL_ARRAY_BUFFER,b.vbo);
    const std::size_t bytes=data.size()*sizeof(glm::vec3);
    if(bytes>b.capacity){b.capacity=std::max<std::size_t>(bytes,4096);glBufferData(GL_ARRAY_BUFFER,static_cast<GLsizeiptr>(b.capacity),nullptr,usage);}
    if(!data.empty()) glBufferSubData(GL_ARRAY_BUFFER,0,static_cast<GLsizeiptr>(bytes),data.data());
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(glm::vec3),nullptr);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

static void DrawPoints(Renderer& r,const std::vector<glm::vec3>& points,float size,const glm::vec3& color,float alpha,const glm::mat4& view,const glm::mat4& projection)
{
    if(points.empty()) return;
    UploadBuffer(r.dynamicPoints,points,GL_STREAM_DRAW);
    glUseProgram(r.pointProgram);
    glUniformMatrix4fv(glGetUniformLocation(r.pointProgram,"uView"),1,GL_FALSE,glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(r.pointProgram,"uProj"),1,GL_FALSE,glm::value_ptr(projection));
    glUniform1f(glGetUniformLocation(r.pointProgram,"uSize"),size);
    glUniform3fv(glGetUniformLocation(r.pointProgram,"uColor"),1,glm::value_ptr(color));
    glUniform1f(glGetUniformLocation(r.pointProgram,"uAlpha"),alpha);
    glBindVertexArray(r.dynamicPoints.vao);
    glDrawArrays(GL_POINTS,0,static_cast<GLsizei>(points.size()));
    glBindVertexArray(0);
}

static void DrawLineStrip(Renderer& r,const std::vector<glm::vec3>& points,float width,const glm::vec3& color,float alpha,const glm::mat4& view,const glm::mat4& projection)
{
    if(points.size()<2) return;
    UploadBuffer(r.dynamicLines,points,GL_STREAM_DRAW);
    glUseProgram(r.lineProgram);
    glUniformMatrix4fv(glGetUniformLocation(r.lineProgram,"uView"),1,GL_FALSE,glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(r.lineProgram,"uProj"),1,GL_FALSE,glm::value_ptr(projection));
    glUniform3fv(glGetUniformLocation(r.lineProgram,"uColor"),1,glm::value_ptr(color));
    glUniform1f(glGetUniformLocation(r.lineProgram,"uAlpha"),alpha);
    glLineWidth(width);
    glBindVertexArray(r.dynamicLines.vao);
    glDrawArrays(GL_LINE_STRIP,0,static_cast<GLsizei>(points.size()));
    glBindVertexArray(0);
}

static void DrawGlowSphere(Renderer& r,const glm::vec3& pos,float scale,const glm::vec3& color,const glm::mat4& view,const glm::mat4& projection,double time,int body)
{
    for(int i=5;i>=1;--i)
    {
        const float s=scale*(1.0f+0.16f*static_cast<float>(i));
        const float emission=0.03f*static_cast<float>(6-i);
        DrawMesh(r,r.sphere,pos,glm::vec3(s),color,emission,0.0f,0.0f,body,view,projection,time);
    }
}

static std::vector<glm::vec3> BuildHaloCurve(const Renderer& r)
{
    std::vector<glm::vec3> curve;
    curve.reserve(361);
    const glm::dvec3 l1{r.points[0].x, r.points[0].y, 0.0};
    for(int i=0;i<=360;++i)
    {
        const double phase = 2.0 * Physics::PI * static_cast<double>(i) / 360.0;
        const glm::dvec3 q = l1 + glm::dvec3(
            0.0016 * std::cos(phase),
            0.0042 * std::sin(phase),
            0.0028 * std::sin(2.0 * phase)
        );
        curve.push_back(InertialDisplayPoint(q, r.simDays));
    }
    return curve;
}

static void BuildSolarWind(const Renderer& r,double time,std::vector<glm::vec3>& out)
{
    out.clear();
    if(!r.showSolarWind) return;

    const glm::vec3 sun=CurrentSunDisplay(r.simDays);
    const int count=static_cast<int>(900.0*r.windDensity);
    out.reserve(static_cast<std::size_t>(count));

    // Full-360-degree 3D radial solar-wind visualization.
    // This is a visualization model, not an MHD plasma solution.
    for(int i=0;i<count;++i)
    {
        const float f=static_cast<float>(i)/static_cast<float>(std::max(1,count));
        const float az=static_cast<float>(2.0*Physics::PI)*std::fmod(f*137.50776f,1.0f);
        const float z=1.0f-2.0f*std::fmod(f*0.61803399f,1.0f);
        const float rr=std::sqrt(std::max(0.0f,1.0f-z*z));
        const float phase=std::fmod(static_cast<float>(time*0.45)+f*3.0f,1.0f);
        const float radius=0.18f+phase*4.6f;
        const float wav=0.035f*std::sin(10.0f*az+static_cast<float>(time)*2.0f+f*30.0f);

        const glm::vec3 dir{rr*std::cos(az),rr*std::sin(az),z};
        const glm::vec3 tangent{-std::sin(az),std::cos(az),0.0f};

        out.push_back(
            sun
            + dir*(radius+wav)
            + tangent*(0.03f*std::sin(static_cast<float>(time)+f*20.0f))
        );
    }
}

static std::vector<std::vector<glm::vec3>> BuildRadiationRays(const Renderer& r,double time)
{
    std::vector<std::vector<glm::vec3>> rays;
    if(!r.showRadiation) return rays;

    const glm::vec3 sun=CurrentSunDisplay(r.simDays);
    constexpr int azimuthCount=48;
    constexpr int elevationCount=7;

    // Full-sphere radiation visualization. Each ray travels radially
    // outward from the Sun, animated by phase along the ray.
    for(int e=0;e<elevationCount;++e)
    {
        const float v=static_cast<float>(e)/static_cast<float>(elevationCount-1);
        const float z=-0.78f+1.56f*v;
        const float rr=std::sqrt(std::max(0.0f,1.0f-z*z));

        for(int a=0;a<azimuthCount;++a)
        {
            const float az=static_cast<float>(2.0*Physics::PI)*static_cast<float>(a)/static_cast<float>(azimuthCount);
            const glm::vec3 dir{rr*std::cos(az),rr*std::sin(az),z};
            std::vector<glm::vec3> line;
            line.reserve(22);

            const float phaseOffset=std::fmod(static_cast<float>(time*0.12)+static_cast<float>(a)*0.013f,1.0f);
            for(int j=0;j<22;++j)
            {
                const float u=std::fmod(static_cast<float>(j)/21.0f+phaseOffset,1.0f);
                line.push_back(sun+dir*(0.18f+u*4.25f));
            }
            rays.push_back(std::move(line));
        }
    }
    return rays;
}

static std::vector<std::vector<glm::vec3>> BuildFlareCurves(const Renderer& r,double time)
{
    std::vector<std::vector<glm::vec3>> flares;
    if(!r.showFlares) return flares;

    const glm::vec3 sun=CurrentSunDisplay(r.simDays);
    const float strength = std::clamp(
        static_cast<float>(r.flareStrength),
        0.10f,
        1.75f
    );

    // Eight active regions distributed around the full 360 degrees.
    for(int k=0;k<8;++k)
    {
        const float base=static_cast<float>(2.0*Physics::PI)*static_cast<float>(k)/8.0f;
        const float pulse=0.50f+0.50f*std::sin(static_cast<float>(time)*1.7f+k*0.73f);
        const float length=0.28f+0.72f*strength*pulse;

        // Curved prominence / flare loop.
        std::vector<glm::vec3> loop;
        loop.reserve(70);
        for(int i=0;i<=68;++i)
        {
            const float u=static_cast<float>(i)/68.0f;
            const float local = base + 0.48f * std::sin(
                u * static_cast<float>(Physics::PI) * 1.5f
                + static_cast<float>(time) * 1.2f
            );

            const float radial =
                0.20f
                + 0.22f * std::sin(u * static_cast<float>(Physics::PI))
                + 0.06f * pulse;
            loop.push_back(sun+glm::vec3(
                radial*std::cos(local),
                radial*std::sin(local),
                0.10f*std::sin(u*Physics::PI*2.0f+k)
            ));
        }
        flares.push_back(std::move(loop));

        // CME / flare jet.
        std::vector<glm::vec3> jet;
        jet.reserve(50);
        glm::vec3 dir{std::cos(base),std::sin(base),0.35f*std::sin(base*2.0f)};
        dir=glm::normalize(dir);
        glm::vec3 side{-dir.y,dir.x,0.0f};
        for(int i=0;i<=48;++i)
        {
            const float u=static_cast<float>(i)/48.0f;
            const float wave=0.055f*std::sin(u*13.0f+static_cast<float>(time)*3.0f+k);
            jet.push_back(sun+dir*(0.12f+length*u)+side*wave*u);
        }
        flares.push_back(std::move(jet));
    }
    return flares;
}

static std::vector<std::vector<glm::vec3>> BuildEarthMagnetic(const Renderer& r,double time)
{
    std::vector<std::vector<glm::vec3>> loops;
    if(!r.showMagnetic) return loops;
    const glm::vec3 earth=CurrentEarthDisplay(r.simDays);
    for(int k=0;k<12;++k)
    {
        const float q=-1.0f+2.0f*static_cast<float>(k)/11.0f;
        std::vector<glm::vec3> line;
        for(int i=0;i<=180;++i)
        {
            const float t=static_cast<float>(i)/180.0f;
            const float a=-0.92f*static_cast<float>(Physics::PI)+1.84f*static_cast<float>(Physics::PI)*t;
            const float r0=0.18f+0.07f*(1.0f-q*q);
            const float nightside=0.70f*(1.0f-std::cos(a))*0.5f;
            line.push_back(earth+glm::vec3(
                r0*std::cos(a)+nightside,
                0.58f*r0*std::sin(a)+0.03f*q,
                0.42f*r0*std::sin(a)*std::sin(a*0.5f)+0.015f*std::sin(time*0.5)
            ));
        }
        loops.push_back(std::move(line));
    }
    return loops;
}

static std::vector<std::vector<glm::vec3>> BuildSolarMagnetic(const Renderer& r,double time)
{
    std::vector<std::vector<glm::vec3>> lines;
    if(!r.showSolarMagnetic) return lines;

    const glm::vec3 sun=CurrentSunDisplay(r.simDays);

    // Full 360-degree 3D dipole-like magnetic visualization. Multiple
    // meridian planes make the field visible from every camera azimuth.
    constexpr int planeCount=20;
    for(int p=0;p<planeCount;++p)
    {
        const float plane=static_cast<float>(2.0*Physics::PI)*static_cast<float>(p)/static_cast<float>(planeCount);
        const float cp=std::cos(plane);
        const float sp=std::sin(plane);

        for(int band=0;band<4;++band)
        {
            const float scale=0.26f+0.06f*static_cast<float>(band);
            std::vector<glm::vec3> line;
            line.reserve(110);

            for(int i=0;i<=108;++i)
            {
                const float u=static_cast<float>(i)/108.0f;
                const float theta=-0.95f*static_cast<float>(Physics::PI)+1.90f*static_cast<float>(Physics::PI)*u;
                const float radial=scale*(1.0f+0.18f*std::cos(theta));
                const float twist=0.10f*std::sin(theta*2.0f+static_cast<float>(time)*0.6f+p);

                const float x=radial*std::cos(theta);
                const float y=radial*std::sin(theta);

                line.push_back(sun+glm::vec3(
                    x,
                    y*cp-twist*sp,
                    y*sp+twist*cp
                ));
            }
            lines.push_back(std::move(line));
        }
    }

    return lines;
}

static glm::vec3 PointMarkerPosition(const Renderer& r,int i)
{
    // The marker itself is ALWAYS the exact calculated point. Labels are
    // offset separately in screen space; the physics marker is never moved.
    return InertialDisplayPoint(NormalizedPoint(r.points[i]), r.simDays);
}

static glm::vec3 PointLeaderTarget(const Renderer& r,int i)
{
    return PointMarkerPosition(r,i);
}

static void DrawLagrangeMarkers(Renderer& r,const glm::mat4& view,const glm::mat4& projection)
{
    if(!r.showPoints) return;

    const std::array<glm::vec3,5> colors={
        glm::vec3(1.0f,0.92f,0.12f),
        glm::vec3(1.0f,0.38f,0.08f),
        glm::vec3(0.45f,0.80f,1.0f),
        glm::vec3(1.0f,0.20f,0.88f),
        glm::vec3(0.25f,1.0f,0.65f)
    };

    const glm::vec3 sun=CurrentSunDisplay(r.simDays);
    const glm::vec3 earth=CurrentEarthDisplay(r.simDays);

    // Geometry guide: all five points are calculated in the rotating
    // CR3BP frame, then transformed to the current inertial display frame.
    glDisable(GL_DEPTH_TEST);

    std::vector<glm::vec3> guide;
    guide.reserve(2);

    // Sun -> Earth baseline.
    guide={sun,earth};
    DrawLineStrip(r,guide,1.0f,{0.55f,0.62f,0.72f},0.28f,view,projection);

    // Sun-Earth-L4-L5 triangle edges.
    const glm::vec3 l4=PointMarkerPosition(r,3);
    const glm::vec3 l5=PointMarkerPosition(r,4);

    guide={sun,l4,earth,l5,sun};
    DrawLineStrip(r,guide,0.9f,{0.46f,0.55f,0.70f},0.24f,view,projection);

    // Exact CR3BP markers. These points are NOT moved for readability.
    for(int i=0;i<5;++i)
    {
        const glm::vec3 actual=PointMarkerPosition(r,i);
        DrawGlowSphere(r,actual,0.020f,colors[i],view,projection,0.0,0);
    }

    // Readability callouts for L1/L2 because a visually enlarged Earth
    // cannot be rendered to physical scale and still remain visible.
    const float calloutX=0.30f;
    const glm::vec3 l1=PointMarkerPosition(r,0);
    const glm::vec3 l2=PointMarkerPosition(r,1);

    const glm::vec3 l1Callout=earth+glm::normalize(l1-earth)*calloutX+glm::vec3(0.0f,0.10f,0.0f);
    const glm::vec3 l2Callout=earth+glm::normalize(l2-earth)*calloutX+glm::vec3(0.0f,-0.10f,0.0f);

    guide={l1,l1Callout};
    DrawLineStrip(r,guide,1.6f,colors[0],0.85f,view,projection);
    guide={l2,l2Callout};
    DrawLineStrip(r,guide,1.6f,colors[1],0.85f,view,projection);

    DrawGlowSphere(r,l1Callout,0.030f,colors[0],view,projection,0.0,0);
    DrawGlowSphere(r,l2Callout,0.030f,colors[1],view,projection,0.0,0);

    glEnable(GL_DEPTH_TEST);
}

static bool WorldToScreen(const glm::vec3& p,const glm::mat4& view,const glm::mat4& proj,int width,int height,ImVec2& out)
{
    const glm::vec4 clip=proj*view*glm::vec4(p,1.0f);
    if(clip.w<=0.0f) return false;
    const glm::vec3 ndc=glm::vec3(clip)/clip.w;
    if(ndc.x<-1.2f||ndc.x>1.2f||ndc.y<-1.2f||ndc.y>1.2f) return false;
    out=ImVec2((ndc.x*0.5f+0.5f)*static_cast<float>(width),(1.0f-(ndc.y*0.5f+0.5f))*static_cast<float>(height));
    return true;
}

static void DrawWorldLabels(const Renderer& r,const glm::mat4& view,const glm::mat4& projection,int width,int height)
{
    if(!r.showWorldLabels||!r.showPoints) return;

    auto* dl=ImGui::GetForegroundDrawList();
    const ImU32 text=IM_COL32(235,245,255,255);
    const std::array<ImU32,5> colors={
        IM_COL32(255,235,50,255),
        IM_COL32(255,120,30,255),
        IM_COL32(130,210,255,255),
        IM_COL32(255,90,220,255),
        IM_COL32(70,255,170,255)
    };

    const glm::vec3 earth=CurrentEarthDisplay(r.simDays);
    const glm::vec3 l1=PointMarkerPosition(r,0);
    const glm::vec3 l2=PointMarkerPosition(r,1);

    const glm::vec3 l1Callout=earth+glm::normalize(l1-earth)*0.30f+glm::vec3(0.0f,0.10f,0.0f);
    const glm::vec3 l2Callout=earth+glm::normalize(l2-earth)*0.30f+glm::vec3(0.0f,-0.10f,0.0f);

    for(int i=0;i<5;++i)
    {
        glm::vec3 labelWorld=PointMarkerPosition(r,i);
        if(i==0) labelWorld=l1Callout;
        if(i==1) labelWorld=l2Callout;

        ImVec2 s;
        if(!WorldToScreen(labelWorld,view,projection,width,height,s)) continue;

        dl->AddCircleFilled(s,5.5f,colors[i]);
        dl->AddCircle(s,8.0f,IM_COL32(255,255,255,180),24,1.0f);
        dl->AddText(ImVec2(s.x+10.0f,s.y-9.0f),text,r.points[i].name.c_str());
    }
}

static void SetupImGui(GLFWwindow* window)
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io=ImGui::GetIO();
    io.ConfigFlags|=ImGuiConfigFlags_NavEnableKeyboard;
    ImGui::StyleColorsDark();
    ImGuiStyle& style=ImGui::GetStyle();
    style.WindowRounding=10.0f;
    style.ChildRounding=7.0f;
    style.FrameRounding=5.0f;
    style.WindowBorderSize=1.0f;
    style.Colors[ImGuiCol_WindowBg]=ImVec4(0.008f,0.012f,0.028f,0.94f);
    style.Colors[ImGuiCol_Header]=ImVec4(0.08f,0.20f,0.42f,0.80f);
    style.Colors[ImGuiCol_HeaderHovered]=ImVec4(0.12f,0.30f,0.58f,0.90f);
    style.Colors[ImGuiCol_CheckMark]=ImVec4(0.30f,0.92f,1.0f,1.0f);
    ImGui_ImplGlfw_InitForOpenGL(window,true);
    ImGui_ImplOpenGL3_Init("#version 330");
}

static void ScrollCallback(GLFWwindow* window,double,double yOffset)
{
    ImGui_ImplGlfw_ScrollCallback(window,0.0,yOffset);
    Renderer* r=static_cast<Renderer*>(glfwGetWindowUserPointer(window));
    if(!r||ImGui::GetIO().WantCaptureMouse) return;
    r->distance*=std::exp(static_cast<float>(-yOffset*0.14));
    r->distance=std::clamp(r->distance,2.0f,20.0f);
}

static void UpdateCameraFromMouse(Renderer& r)
{
    if(ImGui::GetIO().WantCaptureMouse){r.rotatingCamera=false;return;}
    if(glfwGetMouseButton(r.window,GLFW_MOUSE_BUTTON_RIGHT)==GLFW_PRESS)
    {
        double x=0.0,y=0.0;
        glfwGetCursorPos(r.window,&x,&y);
        if(!r.rotatingCamera){r.rotatingCamera=true;r.lastMouseX=x;r.lastMouseY=y;}
        else
        {
            const double dx=x-r.lastMouseX;
            const double dy=y-r.lastMouseY;
            r.yaw+=static_cast<float>(dx*0.006);
            r.pitch+=static_cast<float>(dy*0.006);
            r.pitch=std::clamp(r.pitch,-1.45f,1.45f);
            r.overview=false;
            r.lastMouseX=x;r.lastMouseY=y;
        }
    }
    else r.rotatingCamera=false;
}

static void BuildUI(Renderer& r)
{
    if(!r.showInfo) return;
    ImGui::SetNextWindowSize(ImVec2(455,840),ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(ImVec2(15,15),ImGuiCond_FirstUseEver);
    ImGui::Begin("ADITYA-L1 // FINAL MISSION DIGITAL TWIN",nullptr,ImGuiWindowFlags_NoCollapse);
    ImGui::Text("SUN-EARTH CR3BP | SOLAR OBSERVATORY");
    ImGui::Separator();
    ImGui::Text("Mission: 3-D halo orbit around Sun-Earth L1");
    ImGui::Text("Reference Earth-L1: %.6f M km",r.points[0].earthKm/1.0e6);
    ImGui::Text("Halo period: %.2f days",Physics::ADITYA_HALO_PERIOD_DAYS);
    ImGui::Text("Nominal mission life: %.0f years",Physics::MISSION_LIFETIME_YEARS);
    ImGui::Text("Launch: 02 Sep 2023 | PSLV-C57");
    ImGui::Text("Halo insertion: 06 Jan 2024");
    ImGui::Text("1st halo completion: 02 Jul 2024");

    ImGui::Separator();
    ImGui::Text("SIMULATION");
    ImGui::Checkbox("Pause",&r.paused);
    float speed=static_cast<float>(r.simDaysPerSecond);
    if(ImGui::SliderFloat("Sim days / real second",&speed,0.01f,0.50f,"%.2f")) r.simDaysPerSecond=speed;
    if(ImGui::Button("Reset camera / overview")){r.overview=true;r.yaw=0.0f;r.pitch=0.0f;r.distance=8.5f;}
    ImGui::SameLine();
    ImGui::Checkbox("Follow spacecraft",&r.followSpacecraft);
    ImGui::Checkbox("Show Lagrange points",&r.showPoints);
    ImGui::Checkbox("Show point labels",&r.showWorldLabels);
    ImGui::Checkbox("Show illustrative halo",&r.showHalo);
    ImGui::Checkbox("Show experimental CR3BP path",&r.showTrajectory);
    ImGui::Checkbox("Overview camera",&r.overview);

    ImGui::Separator();
    ImGui::Text("SPACE-WEATHER VISUALIZATIONS");
    ImGui::Checkbox("Solar flares / CME",&r.showFlares);
    ImGui::Checkbox("Solar radiation",&r.showRadiation);
    ImGui::Checkbox("Solar wind",&r.showSolarWind);
    ImGui::Checkbox("Sun magnetic field",&r.showSolarMagnetic);
    ImGui::Checkbox("Earth magnetic field",&r.showMagnetic);
    float fs=static_cast<float>(r.flareStrength);
    if(ImGui::SliderFloat("Flare strength",&fs,0.10f,1.75f,"%.2f")) r.flareStrength=fs;
    float wd=static_cast<float>(r.windDensity);
    if(ImGui::SliderFloat("Solar-wind density",&wd,0.20f,2.50f,"%.2f")) r.windDensity=wd;

    ImGui::Separator();
    ImGui::Text("SPACECRAFT MODE");
    const char* modes[]={"MISSION MODE — L1 HALO","EXPERIMENT — L1","EXPERIMENT — L2","EXPERIMENT — L3","EXPERIMENT — L4","EXPERIMENT — L5","EXPERIMENT — CUSTOM"};
    int combo=r.experimental?r.positionMode+1:0;
    if(ImGui::Combo("Position",&combo,modes,7))
    {
        if(combo==0) r.experimental=false;
        else
        {
            r.experimental=true;
            r.positionMode=combo-1;
            if(r.positionMode>=0&&r.positionMode<5) r.experimentalPosition=NormalizedPoint(r.points[r.positionMode]);
        }
    }
    if(r.experimental&&r.positionMode==5)
    {
        float x=static_cast<float>(r.experimentalPosition.x);
        float y=static_cast<float>(r.experimentalPosition.y);
        float z=static_cast<float>(r.experimentalPosition.z);
        ImGui::DragFloat("Custom X (normalized)",&x,0.0005f,-2.0f,2.0f);
        ImGui::DragFloat("Custom Y (normalized)",&y,0.0005f,-2.0f,2.0f);
        ImGui::DragFloat("Custom Z (normalized)",&z,0.0005f,-1.0f,1.0f);
        r.experimentalPosition={x,y,z};
    }

    ImGui::Separator();
    ImGui::Text("ALL FIVE LAGRANGE POINTS");
    if(ImGui::BeginTable("Lpoints",5,ImGuiTableFlags_Borders|ImGuiTableFlags_RowBg|ImGuiTableFlags_SizingFixedFit))
    {
        ImGui::TableSetupColumn("Point"); ImGui::TableSetupColumn("x"); ImGui::TableSetupColumn("y"); ImGui::TableSetupColumn("Earth km"); ImGui::TableSetupColumn("Sun km"); ImGui::TableHeadersRow();
        for(const auto& p:r.points)
        {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0); ImGui::TextUnformatted(p.name.c_str());
            ImGui::TableSetColumnIndex(1); ImGui::Text("%.9f",p.x);
            ImGui::TableSetColumnIndex(2); ImGui::Text("%.9f",p.y);
            ImGui::TableSetColumnIndex(3); ImGui::Text("%.0f",p.earthKm);
            ImGui::TableSetColumnIndex(4); ImGui::Text("%.0f",p.sunKm);
        }
        ImGui::EndTable();
    }
    static const char* selectedPointLabels[] = {"L1","L2","L3","L4","L5"};
    ImGui::Combo("Selected point",&r.selectedPoint,selectedPointLabels,5);
    const auto& p=r.points[r.selectedPoint];
    ImGui::Text("Selected: %s",p.name.c_str());
    ImGui::Text("Normalized: (%.15f, %.15f, 0)",p.x,p.y);
    ImGui::Text("Earth distance: %,.3f km",p.earthKm);
    ImGui::Text("Sun distance: %,.3f km",p.sunKm);
    ImGui::Text("Stability: %s",r.selectedPoint<3?"dynamically unstable":"linearly stable in CR3BP");
    const glm::dvec3 inertialSelected = RotateNormalized(
        glm::dvec3{p.x,p.y,0.0},
        OrbitalAngleRad(r.simDays)
    );
    ImGui::Text("Current inertial normalized: (%.9f, %.9f, 0)", inertialSelected.x, inertialSelected.y);
    ImGui::TextWrapped("L1-L5 coordinates above are exact for the circular CR3BP rotating frame. The markers are exact; only the L1/L2 callout markers are offset for readability because the Earth model is intentionally enlarged.");

    ImGui::Separator();
    ImGui::Text("LIVE SPACECRAFT");
    const glm::vec3 sc=MissionSpacecraftDisplay(r);
    glm::dvec3 norm;
    if(r.experimental)
        norm = RotateNormalized(r.experimentalPosition, OrbitalAngleRad(r.simDays));
    else
    {
        const double phase = 2.0 * Physics::PI * (r.simDays / Physics::ADITYA_HALO_PERIOD_DAYS);
        norm = glm::dvec3(
            r.points[0].x + 0.0016 * std::cos(phase),
            r.points[0].y + 0.0042 * std::sin(phase),
            0.0028 * std::sin(2.0 * phase)
        );
    }
    ImGui::Text("Simulation day: %.3f",r.simDays);
    ImGui::Text("Earth-Sun orbital angle: %.3f deg",OrbitalAngleRad(r.simDays)*180.0/Physics::PI);
    ImGui::Text("Spacecraft distance from Earth: %,.0f km",EarthDistanceKm(norm));
    ImGui::Text("Spacecraft distance from Sun: %,.0f km",SunDistanceKm(norm));
    ImGui::Text("Solar irradiance estimate: %.2f W/m^2",1361.0*std::pow(Physics::AU_KM/std::max(1.0,SunDistanceKm(norm)),2.0));
    ImGui::Text("Spacecraft reference region: L1");

    if(ImGui::TreeNode("PHYSICS CONSTANTS"))
    {
        ImGui::Text("mu = %.15e",Mu());
        ImGui::Text("1 AU = %.3f km",Physics::AU_KM);
        ImGui::Text("Mean motion = %.12e rad/s",MeanMotion());
        const glm::dvec3 earthRotating{1.0 - Mu(),0.0,0.0};
        const glm::dvec3 earthInertial = RotateNormalized(earthRotating, OrbitalAngleRad(r.simDays));
        ImGui::Text("Earth inertial normalized: (%.9f, %.9f, 0)", earthInertial.x, earthInertial.y);
        ImGui::Text("Velocity scale = %.6f km/s",VelocityScaleKmS());
        ImGui::Text("Time scale = %.6f days",TimeScaleDays());
        ImGui::Text("L1-L5 are exact in the rotating CR3BP frame");
        ImGui::Text("Inertial view rotates them with Earth-Sun mean motion");
        ImGui::Text("Earth/Sun positions therefore change continuously with simulation time");
        ImGui::Text("Jacobi integral conserved by RK4 test path");
        ImGui::TreePop();
    }

    if(ImGui::TreeNode("ADITYA-L1 PAYLOADS — 7"))
    {
        ImGui::BulletText("VELC — Visible Emission Line Coronagraph");
        ImGui::BulletText("SUIT — Solar Ultraviolet Imaging Telescope");
        ImGui::BulletText("SoLEXS — Solar Low Energy X-ray Spectrometer");
        ImGui::BulletText("HEL1OS — High Energy L1 Orbiting X-ray Spectrometer");
        ImGui::BulletText("ASPEX — Aditya Solar wind Particle Experiment");
        ImGui::BulletText("PAPA — Plasma Analyser Package For Aditya");
        ImGui::BulletText("MAG — tri-axial in-situ magnetometer");
        ImGui::TreePop();
    }

    if(ImGui::TreeNode("MISSION / SCIENCE"))
    {
        ImGui::BulletText("Study photosphere, chromosphere and corona");
        ImGui::BulletText("Observe solar flares and CMEs");
        ImGui::BulletText("Measure solar-wind particles and fields before Earth");
        ImGui::BulletText("Continuous solar view from the L1 vantage point");
        ImGui::BulletText("Official halo period reference: 177.86 days");
        ImGui::BulletText("Illustrative graphics are not the operational ephemeris");
        ImGui::BulletText("Solar wind / radiation / flare / field layers are visualization models");
        ImGui::TreePop();
    }

    ImGui::Separator();
    ImGui::TextWrapped("Controls: RMB drag = camera orbit | wheel = zoom | P = pause | I = panel | F = fullscreen | R = reset overview | O = overview camera");
    ImGui::End();
}

static void DrawSystemStatus(Renderer& r,int width)
{
    ImGui::SetNextWindowPos(ImVec2(static_cast<float>(width-350),15),ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(335,215),ImGuiCond_Always);
    ImGui::Begin("SYSTEM STATUS",nullptr,ImGuiWindowFlags_NoResize|ImGuiWindowFlags_NoCollapse);
    ImGui::Text("SUN-EARTH CR3BP");
    ImGui::Text("mu: %.12e",Mu());
    ImGui::Text("L1: %,.0f km Earthward",r.points[0].earthKm);
    ImGui::Text("L2: %,.0f km anti-solar",r.points[1].earthKm);
    ImGui::Text("L3: %,.0f km from Earth",r.points[2].earthKm);
    ImGui::Text("Halo period: %.2f d",Physics::ADITYA_HALO_PERIOD_DAYS);
    ImGui::Text("Simulation: %.2f d",r.simDays);
    ImGui::Text("Camera: %s",r.overview?"OVERVIEW":"3-D FREE ORBIT");
    ImGui::Separator();
    ImGui::Text("Flare %s | Radiation %s",r.showFlares?"ON":"OFF",r.showRadiation?"ON":"OFF");
    ImGui::Text("Wind %s | Earth B %s",r.showSolarWind?"ON":"OFF",r.showMagnetic?"ON":"OFF");
    ImGui::Text("Sun B %s",r.showSolarMagnetic?"ON":"OFF");
    ImGui::End();
}

int main()
{
    try
    {
        std::cout << std::setprecision(15);
        const auto points=ComputeLagrangePoints();

        std::cout << "ADITYA-L1 // FINAL MISSION DIGITAL TWIN\n";
        std::cout << "mu = " << Mu() << "\n";
        for(const auto& p:points)
            std::cout << p.name << " = (" << p.x << ", " << p.y << ") | Earth = " << p.earthKm << " km | Sun = " << p.sunKm << " km\n";

        if(!glfwInit()) return 1;
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
        glfwWindowHint(GLFW_OPENGL_PROFILE,GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_SAMPLES,4);

        GLFWwindow* window=glfwCreateWindow(1680,980,"Aditya-L1 // Final Mission Digital Twin",nullptr,nullptr);
        if(!window){glfwTerminate();return 1;}
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1);

        if(gladLoadGL(glfwGetProcAddress)==0){glfwDestroyWindow(window);glfwTerminate();return 1;}

        Renderer r;
        r.window=window;
        glfwSetWindowUserPointer(window,&r);
        r.points=points;
        r.experimentalPosition=NormalizedPoint(r.points[0]);

        r.meshProgram=CreateProgram(MeshVS,MeshFS);
        r.pointProgram=CreateProgram(PointVS,PointFS);
        r.lineProgram=CreateProgram(LineVS,LineFS);
        if(!r.meshProgram||!r.pointProgram||!r.lineProgram){glfwDestroyWindow(window);glfwTerminate();return 1;}

        r.sphere=MakeSphere(72,144);
        r.spacecraftBody=MakeCube();
        r.spacecraftPanel=MakeCube();

        std::mt19937 rng(20260813u);
        std::uniform_real_distribution<float> sx(-18.0f,18.0f), sy(-10.0f,10.0f), sz(-18.0f,18.0f);
        r.stars.reserve(7000);
        for(int i=0;i<7000;++i) r.stars.push_back({sx(rng),sy(rng),sz(rng)});

        SetupImGui(window);
        glfwSetScrollCallback(window,ScrollCallback);

        glEnable(GL_DEPTH_TEST);
        glEnable(GL_PROGRAM_POINT_SIZE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_MULTISAMPLE);
        glClearColor(0.0015f,0.003f,0.010f,1.0f);

        bool prevP=false,prevI=false,prevF11=false,prevR=false,prevO=false;
        double lastTime=glfwGetTime();

        while(!glfwWindowShouldClose(window))
        {
            glfwPollEvents();
            const double now=glfwGetTime();
            const double dt=std::min(0.1,now-lastTime);
            lastTime=now;

            const bool pKey=glfwGetKey(window,GLFW_KEY_P)==GLFW_PRESS;
            const bool iKey=glfwGetKey(window,GLFW_KEY_I)==GLFW_PRESS;
            const bool f11Key=glfwGetKey(window,GLFW_KEY_F11)==GLFW_PRESS;
            const bool rKey=glfwGetKey(window,GLFW_KEY_R)==GLFW_PRESS;
            const bool oKey=glfwGetKey(window,GLFW_KEY_O)==GLFW_PRESS;

            if(!ImGui::GetIO().WantCaptureKeyboard)
            {
                if(pKey&&!prevP) r.paused=!r.paused;
                if(iKey&&!prevI) r.showInfo=!r.showInfo;
                if(rKey&&!prevR){r.overview=true;r.yaw=0.0f;r.pitch=0.0f;r.distance=8.5f;r.followSpacecraft=false;}
                if(oKey&&!prevO) r.overview=!r.overview;
            }

            if(f11Key&&!prevF11)
            {
                const bool fullscreen=glfwGetWindowMonitor(window)!=nullptr;
                if(!fullscreen)
                {
                    GLFWmonitor* monitor=glfwGetPrimaryMonitor();
                    const GLFWvidmode* mode=glfwGetVideoMode(monitor);
                    glfwSetWindowMonitor(window,monitor,0,0,mode->width,mode->height,mode->refreshRate);
                }
                else glfwSetWindowMonitor(window,nullptr,100,80,1680,980,0);
            }

            prevP=pKey;prevI=iKey;prevF11=f11Key;prevR=rKey;prevO=oKey;

            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            UpdateCameraFromMouse(r);
            if(!r.paused)
            {
                r.simDays+=r.simDaysPerSecond*dt;
                if(r.simDays>=Physics::ADITYA_HALO_PERIOD_DAYS) r.simDays=std::fmod(r.simDays,Physics::ADITYA_HALO_PERIOD_DAYS);
            }

            int width=0,height=0;
            glfwGetFramebufferSize(window,&width,&height);
            glViewport(0,0,width,height);
            glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

            const glm::vec3 spacecraft=MissionSpacecraftDisplay(r);
            glm::vec3 target(0.0f);
            if(r.followSpacecraft) target=spacecraft;

            glm::mat4 view;
            glm::mat4 projection;

            if(r.overview)
            {
                const glm::vec3 camera(0.0f,0.0f,r.distance);
                view=glm::lookAt(camera,target,glm::vec3(0.0f,1.0f,0.0f));
                const float aspect=static_cast<float>(width)/static_cast<float>(std::max(1,height));
                const float viewHeight=7.4f*std::clamp(r.distance/8.5f,0.65f,1.65f);
                projection=glm::ortho(-viewHeight*aspect*0.5f,viewHeight*aspect*0.5f,-viewHeight*0.5f,viewHeight*0.5f,-40.0f,40.0f);
            }
            else
            {
                const glm::vec3 camera=CameraPosition(r);
                view=glm::lookAt(camera+target,target,glm::vec3(0.0f,1.0f,0.0f));
                const float aspect=static_cast<float>(width)/static_cast<float>(std::max(1,height));
                projection=glm::perspective(glm::radians(48.0f),aspect,0.01f,100.0f);
            }

            // Stars
            DrawPoints(r,r.stars,1.6f,{0.70f,0.82f,1.0f},0.70f,view,projection);

            // Time-dependent inertial Sun/Earth positions from the exact
            // circular Sun-Earth CR3BP rotation.
            const glm::vec3 sunPos=CurrentSunDisplay(r.simDays);
            const glm::vec3 earthPos=CurrentEarthDisplay(r.simDays);

            // Sun
            DrawGlowSphere(r,sunPos,0.42f,{1.0f,0.25f,0.01f},view,projection,now,1);
            DrawMesh(r,r.sphere,sunPos,{0.31f,0.31f,0.31f},{1.0f,0.43f,0.03f},0.70f,0.18f,static_cast<float>(now*0.04),1,view,projection,now);

            // Earth
            DrawGlowSphere(r,earthPos,0.14f,{0.03f,0.15f,0.90f},view,projection,now,2);
            DrawMesh(r,r.sphere,earthPos,{0.105f,0.105f,0.105f},{0.03f,0.18f,0.60f},0.02f,0.35f,static_cast<float>(now*0.12),2,view,projection,now);

            // Lagrange points: exact physics positions + readability markers with leaders.
            DrawLagrangeMarkers(r,view,projection);

            // Halo envelope
            if(r.showHalo) DrawLineStrip(r,BuildHaloCurve(r),2.0f,{0.08f,0.84f,1.0f},0.80f,view,projection);

            // Spacecraft: one spacecraft only.
            if(!r.experimental)
            {
                DrawMesh(r,r.spacecraftBody,spacecraft,{0.085f,0.085f,0.085f},{0.86f,0.90f,0.98f},0.22f,0.45f,static_cast<float>(now*0.35),0,view,projection,now);
                DrawMesh(r,r.spacecraftPanel,spacecraft+glm::vec3(-0.16f,0.0f,0.0f),{0.11f,0.018f,0.052f},{0.02f,0.16f,0.95f},0.03f,0.10f,static_cast<float>(now*0.05),0,view,projection,now);
                DrawMesh(r,r.spacecraftPanel,spacecraft+glm::vec3(0.16f,0.0f,0.0f),{0.11f,0.018f,0.052f},{0.02f,0.16f,0.95f},0.03f,0.10f,static_cast<float>(now*0.05),0,view,projection,now);
                DrawGlowSphere(r,spacecraft,0.065f,{0.05f,1.0f,1.0f},view,projection,now,0);
            }
            else
            {
                DrawGlowSphere(r,spacecraft,0.09f,{1.0f,0.20f,0.85f},view,projection,now,0);
                DrawMesh(r,r.spacecraftBody,spacecraft,{0.09f,0.09f,0.09f},{0.96f,0.96f,1.0f},0.18f,0.45f,static_cast<float>(now*0.22),0,view,projection,now);
            }

            if(r.showTrajectory)
            {
                if(r.trajectory.empty())
                {
                    Physics::State s{{r.points[0].x-0.001,0.0,0.001},{0.0,0.002,0.0}};
                    constexpr int N=18000;
                    constexpr double dtSim=0.0010;
                    r.trajectory.reserve(N);
                    for(int i=0;i<N;++i){r.trajectory.push_back(DisplayPoint(s.r));s=RK4Step(s,dtSim);}
                }
                DrawLineStrip(r,r.trajectory,1.3f,{0.08f,0.58f,1.0f},0.28f,view,projection);
            }

            // Solar wind
            std::vector<glm::vec3> wind;
            BuildSolarWind(r,now,wind);
            DrawPoints(r,wind,2.2f,{0.20f,0.78f,1.0f},0.52f,view,projection);

            // Solar radiation
            const auto radiation=BuildRadiationRays(r,now);
            for(const auto& ray:radiation) DrawLineStrip(r,ray,1.0f,{1.0f,0.74f,0.20f},static_cast<float>(0.26*r.radiationStrength),view,projection);

            // Flares / CME
            const auto flares=BuildFlareCurves(r,now);
            for(const auto& flare:flares) DrawLineStrip(r,flare,2.8f,{1.0f,0.18f,0.03f},static_cast<float>(0.70*r.flareStrength),view,projection);

            // Sun magnetic field
            const auto sunB=BuildSolarMagnetic(r,now);
            for(const auto& loop:sunB) DrawLineStrip(r,loop,1.15f,{1.0f,0.33f,0.08f},0.30f,view,projection);

            // Earth magnetic field
            const auto earthB=BuildEarthMagnetic(r,now);
            for(const auto& loop:earthB) DrawLineStrip(r,loop,1.25f,{0.12f,0.52f,1.0f},0.38f,view,projection);

            BuildUI(r);
            DrawSystemStatus(r,width);
            DrawWorldLabels(r,view,projection,width,height);

            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            glfwSwapBuffers(window);
        }

        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        DestroyMesh(r.sphere);
        DestroyMesh(r.spacecraftBody);
        DestroyMesh(r.spacecraftPanel);
        DestroyLineBuffer(r.dynamicLines);
        DestroyLineBuffer(r.dynamicPoints);
        if(r.meshProgram) glDeleteProgram(r.meshProgram);
        if(r.pointProgram) glDeleteProgram(r.pointProgram);
        if(r.lineProgram) glDeleteProgram(r.lineProgram);
        glfwDestroyWindow(window);
        glfwTerminate();
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << '\n';
        return 1;
    }
    catch(...)
    {
        std::cerr << "FATAL: unknown error\n";
        return 1;
    }
}
