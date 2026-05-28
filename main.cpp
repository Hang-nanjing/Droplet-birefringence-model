#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// 定义二维向量结构体和相关函数
struct Vec2 {
    double x{};
    double y{};

    Vec2 operator+(const Vec2& other) const {
        return {x + other.x, y + other.y};
    }

    Vec2 operator-(const Vec2& other) const {
        return {x - other.x, y - other.y};
    }

    Vec2 operator*(double scale) const {
        return {x * scale, y * scale};
    }
};

// 向量点积
double dot(const Vec2& a, const Vec2& b) {
    return a.x * b.x + a.y * b.y;
}

// 向量长度
double norm(const Vec2& v) {
    return std::sqrt(dot(v, v));
}

double clampValue(double value, double low, double high) {
    return std::max(low, std::min(value, high));
}

double radiansToDegrees(double radians) {
    constexpr double pi = 3.14159265358979323846;
    return radians * 180.0 / pi;
}

double directionAngle(const Vec2& direction) {
    return std::atan2(direction.y, direction.x);
}

double shortestAngleSweep(double from, double to) {
    constexpr double pi = 3.14159265358979323846;
    double sweep = to - from;
    while (sweep > pi) {
        sweep -= 2.0 * pi;
    }
    while (sweep < -pi) {
        sweep += 2.0 * pi;
    }
    return sweep;
}

// 向量归一化
Vec2 normalize(const Vec2& v) {
    const double length = norm(v);
    if (length == 0.0) {
        throw std::runtime_error("Cannot normalize a zero-length vector.");
    }
    return {v.x / length, v.y / length};
}

double angleBetweenUnitDirectionsDeg(const Vec2& a, const Vec2& b) {
    const double c = clampValue(dot(normalize(a), normalize(b)), -1.0, 1.0);
    return radiansToDegrees(std::acos(c));
}

// 格式化数字为字符串，处理非有限值
std::string formatNumber(double value) {
    if (!std::isfinite(value)) {
        return "";
    }
    std::ostringstream out;
    out << std::fixed << std::setprecision(9) << value;
    return out.str();
} 

std::string formatAngle(double value) {
    if (!std::isfinite(value)) {
        return "";
    }
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value << " deg";
    return out.str();
}

// 定义系统参数
struct SystemParameters {
    double backgroundX = 0.0;
    double backgroundSizeY = 15.0;

    Vec2 cameraCenter{200.391, 0.0};
    double focalLength = 120.235;

    Vec2 sensorCenter{320.626, 0.0};
    double sensorSizeY = 8.445;

    double pixelSizeMm = 0.00345;
    int sampleStepPixels = 100;

    Vec2 dropletCenter{2.0, 0.0};
    double dropletRadius = 1.5;
    double refractiveIndexAir = 1.0;
    double refractiveIndexWater = 1.33;
};

// 输出文件路径和名称
struct OutputFiles {
    std::string noDropletSvg = "no_droplet_trace.svg";
    std::string dropletSvg = "droplet_trace.svg";
    std::string comparisonSvg = "ray_trace_comparison.svg";
    std::string nearFieldSvg = "near_field_trace.svg";
    std::string displacementCsv = "displacement_table.csv";
};

// svg 可视化设置
struct VisualizationSettings {
    double fullViewHeight = 620.0;
    double nearFieldHeight = 620.0;
    double nearFieldXMin = -1.0;
    double nearFieldXMax = 4.5;
    double nearFieldYMin = -2.2;
    double nearFieldYMax = 2.2;
};

// 定义每条光线追踪的结果结构体
struct TraceResult {
    int index{};
    Vec2 sensorPoint{};
    Vec2 initialDirection{};
    Vec2 q0{};
    Vec2 q{};
    Vec2 delta{};
    Vec2 p1{std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
    Vec2 p2{std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
    Vec2 p1Normal{std::numeric_limits<double>::quiet_NaN(),
                  std::numeric_limits<double>::quiet_NaN()};
    Vec2 p2Normal{std::numeric_limits<double>::quiet_NaN(),
                  std::numeric_limits<double>::quiet_NaN()};
    Vec2 internalDirection{std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN()};
    Vec2 exitDirection{std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::quiet_NaN()};
    double entryIncidentAngleDeg = std::numeric_limits<double>::quiet_NaN();
    double entryRefractedAngleDeg = std::numeric_limits<double>::quiet_NaN();
    double exitIncidentAngleDeg = std::numeric_limits<double>::quiet_NaN();
    double exitRefractedAngleDeg = std::numeric_limits<double>::quiet_NaN();
    std::vector<Vec2> noDropletPath;
    std::vector<Vec2> dropletPath;
    std::string status;
};

// 统计结果摘要
struct TraceSummary {
    int total{};
    int noHit{};
    int refracted{};
    int tir{};
    int invalid{};
    double minDeltaY{std::numeric_limits<double>::infinity()};
    double maxDeltaY{-std::numeric_limits<double>::infinity()};
    double minDeltaMagnitude{std::numeric_limits<double>::infinity()};
    double maxDeltaMagnitude{-std::numeric_limits<double>::infinity()};
};

// 折射函数，计算入射光线经过界面折射后的方向
bool refract(const Vec2& incident,
             const Vec2& normalToIncidentMedium,
             double n1,
             double n2,
             Vec2& transmitted) {
    const Vec2 i = normalize(incident);
    const Vec2 n = normalize(normalToIncidentMedium);
    const double eta = n1 / n2;
    const double ci = -dot(n, i);
    const double k = 1.0 - eta * eta * (1.0 - ci * ci);
    if (k < 0.0) {
        return false;
    }
    transmitted = normalize(i * eta + n * (eta * ci - std::sqrt(k)));
    return true;
}

// 计算光线与竖直背景板的交点
bool intersectRayWithVerticalLine(const Vec2& origin,
                                  const Vec2& direction,
                                  double x,
                                  Vec2& hit) {
    if (std::abs(direction.x) < 1e-12) {
        return false;
    }
    const double t = (x - origin.x) / direction.x;
    if (t < 0.0) {
        return false;
    }
    hit = origin + direction * t;
    return true;
}

// 计算光线与圆形液滴的交点
std::vector<double> intersectRayCircleParameters(const Vec2& origin,
                                                 const Vec2& direction,
                                                 const Vec2& center,
                                                 double radius,
                                                 double minT) {
    const Vec2 oc = origin - center;
    const double a = dot(direction, direction);
    const double b = 2.0 * dot(oc, direction);
    const double c = dot(oc, oc) - radius * radius;
    const double disc = b * b - 4.0 * a * c;

    std::vector<double> hits;
    if (disc < -1e-12) {
        return hits;
    }

    const double clampedDisc = std::max(0.0, disc);
    const double sqrtDisc = std::sqrt(clampedDisc);
    const double invDenom = 1.0 / (2.0 * a);
    const double t0 = (-b - sqrtDisc) * invDenom;
    const double t1 = (-b + sqrtDisc) * invDenom;

    if (t0 > minT) {
        hits.push_back(t0);
    }
    if (t1 > minT && std::abs(t1 - t0) > 1e-9) {
        hits.push_back(t1);
    }
    std::sort(hits.begin(), hits.end());
    return hits;
}

// 生成传感器上的采样点列表
std::vector<Vec2> makeSensorPoints(const SystemParameters& params) {
    const double halfY = params.sensorSizeY * 0.5;
    const double stepMm = params.pixelSizeMm * params.sampleStepPixels;

    std::vector<Vec2> points;
    points.push_back(params.sensorCenter);
    for (double offset = stepMm; offset <= halfY + 1e-9; offset += stepMm) {
        points.push_back({params.sensorCenter.x, params.sensorCenter.y + offset});
        points.push_back({params.sensorCenter.x, params.sensorCenter.y - offset});
    }

    std::sort(points.begin(), points.end(), [](const Vec2& a, const Vec2& b) {
        return a.y > b.y;
    });
    return points;
}

//  对单个传感器点进行光线追踪
TraceResult traceSample(const SystemParameters& params, int index, const Vec2& sensorPoint) {
    TraceResult result;
    result.index = index;
    result.sensorPoint = sensorPoint;
    result.initialDirection = normalize(params.cameraCenter - sensorPoint);

    if (!intersectRayWithVerticalLine(params.cameraCenter,
                                      result.initialDirection,
                                      params.backgroundX,
                                      result.q0)) {
        result.status = "invalid";
        result.q = result.q0;
        result.delta = {0.0, 0.0};
        result.noDropletPath = {sensorPoint, params.cameraCenter};
        result.dropletPath = result.noDropletPath;
        return result;
    }

    result.noDropletPath = {sensorPoint, params.cameraCenter, result.q0};

    const auto entryHits = intersectRayCircleParameters(params.cameraCenter,
                                                       result.initialDirection,
                                                       params.dropletCenter,
                                                       params.dropletRadius,
                                                       1e-9);
    if (entryHits.empty()) {
        result.status = "no_hit";
        result.q = result.q0;
        result.delta = {0.0, 0.0};
        result.dropletPath = result.noDropletPath;
        return result;
    }

    const double tBackground = (params.backgroundX - params.cameraCenter.x) /
                               result.initialDirection.x;
    if (entryHits.front() >= tBackground) {
        result.status = "no_hit";
        result.q = result.q0;
        result.delta = {0.0, 0.0};
        result.dropletPath = result.noDropletPath;
        return result;
    }

    result.p1 = params.cameraCenter + result.initialDirection * entryHits.front();
    const Vec2 n1 = normalize(result.p1 - params.dropletCenter);
    result.p1Normal = n1;

    Vec2 internalDirection;
    if (!refract(result.initialDirection,
                 n1,
                 params.refractiveIndexAir,
                 params.refractiveIndexWater,
                 internalDirection)) {
        result.status = "TIR";
        result.q = result.q0;
        result.delta = {0.0, 0.0};
        result.dropletPath = {sensorPoint, params.cameraCenter, result.p1};
        return result;
    }
    result.internalDirection = internalDirection;
    result.entryIncidentAngleDeg = angleBetweenUnitDirectionsDeg(result.initialDirection * -1.0, n1);
    result.entryRefractedAngleDeg = angleBetweenUnitDirectionsDeg(internalDirection, n1 * -1.0);

    const auto exitHits = intersectRayCircleParameters(result.p1,
                                                      internalDirection,
                                                      params.dropletCenter,
                                                      params.dropletRadius,
                                                      1e-7);
    if (exitHits.empty()) {
        result.status = "invalid";
        result.q = result.q0;
        result.delta = {0.0, 0.0};
        result.dropletPath = {sensorPoint, params.cameraCenter, result.p1};
        return result;
    }

    result.p2 = result.p1 + internalDirection * exitHits.front();
    const Vec2 n2Outward = normalize(result.p2 - params.dropletCenter);
    result.p2Normal = n2Outward;
    const Vec2 n2ToIncidentMedium = n2Outward * -1.0;

    Vec2 exitDirection;
    if (!refract(internalDirection,
                 n2ToIncidentMedium,
                 params.refractiveIndexWater,
                 params.refractiveIndexAir,
                 exitDirection)) {
        result.status = "TIR";
        result.q = result.q0;
        result.delta = {0.0, 0.0};
        result.dropletPath = {sensorPoint, params.cameraCenter, result.p1, result.p2};
        return result;
    }
    result.exitDirection = exitDirection;
    result.exitIncidentAngleDeg = angleBetweenUnitDirectionsDeg(internalDirection * -1.0,
                                                                n2ToIncidentMedium);
    result.exitRefractedAngleDeg = angleBetweenUnitDirectionsDeg(exitDirection, n2Outward);

    if (!intersectRayWithVerticalLine(result.p2,
                                      exitDirection,
                                      params.backgroundX,
                                      result.q)) {
        result.status = "invalid";
        result.q = result.q0;
        result.delta = {0.0, 0.0};
        result.dropletPath = {sensorPoint, params.cameraCenter, result.p1, result.p2};
        return result;
    }

    result.status = "refracted";
    result.delta = result.q - result.q0;
    result.dropletPath = {sensorPoint, params.cameraCenter, result.p1, result.p2, result.q};
    return result;
}

std::vector<std::string> makePalette() {
    return {
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#393b79",
        "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd",
        "#e6550d", "#31a354", "#756bb1", "#636363", "#9ecae1",
        "#fdae6b", "#a1d99b", "#bcbddc", "#969696", "#e41a1c",
    };
}

class SvgWriter {
public:
    SvgWriter(const SystemParameters& params, double height)
        : SvgWriter(params,
                    height,
                    params.backgroundX - 8.0,
                    params.sensorCenter.x + 18.0,
                    -8.5,
                    8.5) {}

    SvgWriter(const SystemParameters& params,
              double height,
              double xMinMm,
              double xMaxMm,
              double yMinMm,
              double yMaxMm)
        : params_(params),
          height_(height),
          xMinMm_(xMinMm),
          xMaxMm_(xMaxMm),
          yMinMm_(yMinMm),
          yMaxMm_(yMaxMm) {
        const double yRangeMm = yMaxMm_ - yMinMm_;
        pixelsPerMm_ = (height_ - 2.0 * margin_) / yRangeMm;
        width_ = (xMaxMm_ - xMinMm_) * pixelsPerMm_ + 2.0 * margin_;
    }

    double width() const {
        return width_;
    }

    double sx(double x) const {
        return margin_ + (x - xMinMm_) * pixelsPerMm_;
    }

    double sy(double y) const {
        return margin_ + (yMaxMm_ - y) * pixelsPerMm_;
    }

    void writeHeader(std::ofstream& svg, const std::string& title) const {
        svg << std::fixed << std::setprecision(3);
        svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width_
            << "\" height=\"" << height_ << "\" viewBox=\"0 0 " << width_
            << " " << height_ << "\">\n";
        svg << "  <rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
        svg << "  <style>"
            << "text{font-family:Arial,Helvetica,sans-serif;font-size:15px;fill:#1f2933;}"
            << ".label{font-size:17px;font-weight:600;}"
            << ".small{font-size:12px;fill:#52606d;}"
            << "</style>\n";
        svg << "  <text class=\"label\" x=\"32\" y=\"34\">" << title << "</text>\n";
    }

    void writeScene(std::ofstream& svg, bool includeDroplet) const {
        const double bgHalf = params_.backgroundSizeY * 0.5;
        const double sensorHalf = params_.sensorSizeY * 0.5;

        svg << "  <line x1=\"" << sx(params_.backgroundX) << "\" y1=\"" << sy(-bgHalf)
            << "\" x2=\"" << sx(params_.backgroundX) << "\" y2=\"" << sy(bgHalf)
            << "\" stroke=\"#111111\" stroke-width=\"5\"/>\n";
        svg << "  <text class=\"small\" x=\"" << sx(params_.backgroundX) + 8.0
            << "\" y=\"" << sy(bgHalf) - 10.0 << "\">background x=0</text>\n";

        svg << "  <line x1=\"" << sx(params_.sensorCenter.x) << "\" y1=\"" << sy(-sensorHalf)
            << "\" x2=\"" << sx(params_.sensorCenter.x) << "\" y2=\"" << sy(sensorHalf)
            << "\" stroke=\"#0f3b82\" stroke-width=\"5\"/>\n";
        svg << "  <text class=\"small\" x=\"" << sx(params_.sensorCenter.x) - 72.0
            << "\" y=\"" << sy(-sensorHalf) + 32.0 << "\">sensor plane</text>\n";

        svg << "  <line x1=\"" << sx(params_.backgroundX - 4.0) << "\" y1=\"" << sy(0.0)
            << "\" x2=\"" << sx(params_.sensorCenter.x + 8.0) << "\" y2=\"" << sy(0.0)
            << "\" stroke=\"#9aa5b1\" stroke-width=\"1.4\"/>\n";

        svg << "  <circle cx=\"" << sx(params_.cameraCenter.x) << "\" cy=\"" << sy(params_.cameraCenter.y)
            << "\" r=\"6\" fill=\"#111111\"/>\n";
        svg << "  <text class=\"small\" x=\"" << sx(params_.cameraCenter.x) - 42.0
            << "\" y=\"" << sy(params_.cameraCenter.y) - 14.0 << "\">camera center</text>\n";

        if (includeDroplet) {
            const double r = params_.dropletRadius * pixelsPerMm_;
            svg << "  <circle cx=\"" << sx(params_.dropletCenter.x) << "\" cy=\""
                << sy(params_.dropletCenter.y) << "\" r=\"" << r
                << "\" fill=\"#dbeafe\" fill-opacity=\"0.45\" stroke=\"#2563eb\""
                << " stroke-width=\"2\"/>\n";
            svg << "  <text class=\"small\" x=\"" << sx(params_.dropletCenter.x) - 42.0
                << "\" y=\"" << sy(params_.dropletCenter.y + params_.dropletRadius) - 8.0
                << "\">droplet R=1.5</text>\n";
        }
    }

    void writePath(std::ofstream& svg,
                   const std::vector<Vec2>& points,
                   const std::string& color,
                   bool dashed,
                   double width) const {
        if (points.size() < 2) {
            return;
        }
        svg << "  <polyline points=\"";
        for (const auto& point : points) {
            svg << sx(point.x) << "," << sy(point.y) << " ";
        }
        svg << "\" fill=\"none\" stroke=\"" << color << "\" stroke-width=\"" << width
            << "\" stroke-linecap=\"round\" stroke-linejoin=\"round\"";
        if (dashed) {
            svg << " stroke-dasharray=\"8 7\"";
        }
        svg << "/>\n";
    }

    void writeSegment(std::ofstream& svg,
                      const Vec2& a,
                      const Vec2& b,
                      const std::string& color,
                      double width) const {
        svg << "  <line x1=\"" << sx(a.x) << "\" y1=\"" << sy(a.y)
            << "\" x2=\"" << sx(b.x) << "\" y2=\"" << sy(b.y)
            << "\" stroke=\"" << color << "\" stroke-width=\"" << width
            << "\" stroke-linecap=\"round\"/>\n";
    }

    void writeText(std::ofstream& svg,
                   const Vec2& position,
                   const std::string& text,
                   const std::string& color,
                   int fontSize) const {
        svg << "  <text x=\"" << sx(position.x) << "\" y=\"" << sy(position.y)
            << "\" font-size=\"" << fontSize << "\" fill=\"" << color
            << "\">" << text << "</text>\n";
    }

    void writeAngleArc(std::ofstream& svg,
                       const Vec2& center,
                       const Vec2& directionA,
                       const Vec2& directionB,
                       double radiusMm,
                       double labelOffsetMm,
                       const std::string& label,
                       const std::string& color,
                       int fontSize) const {
        if (label.empty()) {
            return;
        }

        const double angleA = directionAngle(directionA);
        const double sweep = shortestAngleSweep(angleA, directionAngle(directionB));
        if (std::abs(sweep) < 0.008) {
            return;
        }

        constexpr int segments = 16;
        svg << "  <polyline points=\"";
        for (int i = 0; i <= segments; ++i) {
            const double t = static_cast<double>(i) / segments;
            const double angle = angleA + sweep * t;
            const Vec2 point{
                center.x + radiusMm * std::cos(angle),
                center.y + radiusMm * std::sin(angle),
            };
            svg << sx(point.x) << "," << sy(point.y) << " ";
        }
        svg << "\" fill=\"none\" stroke=\"" << color
            << "\" stroke-width=\"1\" stroke-linecap=\"round\"/>\n";

        const double labelAngle = angleA + sweep * 0.5;
        const double labelRadius = radiusMm + labelOffsetMm;
        const Vec2 labelPosition{
            center.x + labelRadius * std::cos(labelAngle),
            center.y + labelRadius * std::sin(labelAngle),
        };
        writeText(svg, labelPosition, label, color, fontSize);
    }

    void writeFooter(std::ofstream& svg) const {
        svg << "</svg>\n";
    }

private:
    const SystemParameters& params_;
    double width_{};
    double height_{};
    double pixelsPerMm_{};
    double margin_ = 76.0;
    double xMinMm_{};
    double xMaxMm_{};
    double yMinMm_{};
    double yMaxMm_{};
};

void writeTraceSvg(const SystemParameters& params,
                   const VisualizationSettings& visualization,
                   const std::vector<TraceResult>& results,
                   const std::string& fileName,
                   const std::string& title,
                   bool includeDroplet,
                   bool useDropletPaths,
                   bool dashed) {
    std::ofstream svg(fileName);
    if (!svg) {
        throw std::runtime_error("Failed to open SVG output file: " + fileName);
    }

    const auto palette = makePalette();
    SvgWriter writer(params, visualization.fullViewHeight);
    writer.writeHeader(svg, title);
    writer.writeScene(svg, includeDroplet);
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto& path = useDropletPaths ? results[i].dropletPath : results[i].noDropletPath;
        writer.writePath(svg, path, palette[i % palette.size()], dashed, 1.9);
    }
    writer.writeFooter(svg);
}

void writeComparisonSvg(const SystemParameters& params,
                        const VisualizationSettings& visualization,
                        const std::vector<TraceResult>& results,
                        const std::string& fileName) {
    std::ofstream svg(fileName);
    if (!svg) {
        throw std::runtime_error("Failed to open SVG output file: " + fileName);
    }

    const auto palette = makePalette();
    SvgWriter writer(params, visualization.fullViewHeight);
    writer.writeHeader(svg, "Ray Trace Comparison: dashed=no droplet, solid=with droplet");
    writer.writeScene(svg, true);
    for (std::size_t i = 0; i < results.size(); ++i) {
        const std::string& color = palette[i % palette.size()];
        writer.writePath(svg, results[i].noDropletPath, color, true, 1.5);
        writer.writePath(svg, results[i].dropletPath, color, false, 2.0);
    }
    writer.writeFooter(svg);
}

void writeRefractionAnnotations(std::ofstream& svg,
                                const SvgWriter& writer,
                                const TraceResult& result,
                                const std::string& color) {
    if (result.status != "refracted") {
        return;
    }

    const double normalHalfLengthMm = 0.55;
    const double incidentArcRadiusMm = 0.27;
    const double refractedArcRadiusMm = 0.42;
    const int angleFontSize = 6;

    const Vec2 p1Normal = normalize(result.p1Normal);
    const Vec2 p2Normal = normalize(result.p2Normal);

    writer.writeSegment(svg,
                        result.p1 - p1Normal * normalHalfLengthMm,
                        result.p1 + p1Normal * normalHalfLengthMm,
                        color,
                        0.9);
    writer.writeAngleArc(svg,
                         result.p1,
                         result.initialDirection * -1.0,
                         p1Normal,
                         incidentArcRadiusMm,
                         0.30,
                         formatAngle(result.entryIncidentAngleDeg),
                         color,
                         angleFontSize);
    writer.writeAngleArc(svg,
                         result.p1,
                         result.internalDirection,
                         p1Normal * -1.0,
                         refractedArcRadiusMm,
                         0.24,
                         formatAngle(result.entryRefractedAngleDeg),
                         color,
                         angleFontSize);

    writer.writeSegment(svg,
                        result.p2 - p2Normal * normalHalfLengthMm,
                        result.p2 + p2Normal * normalHalfLengthMm,
                        color,
                        0.9);
    writer.writeAngleArc(svg,
                         result.p2,
                         result.internalDirection * -1.0,
                         p2Normal * -1.0,
                         incidentArcRadiusMm,
                         0.38,
                         formatAngle(result.exitIncidentAngleDeg),
                         color,
                         angleFontSize);
    writer.writeAngleArc(svg,
                         result.p2,
                         result.exitDirection,
                         p2Normal,
                         refractedArcRadiusMm,
                         0.34,
                         formatAngle(result.exitRefractedAngleDeg),
                         color,
                         angleFontSize);
}

void writeNearFieldSvg(const SystemParameters& params,
                       const VisualizationSettings& visualization,
                       const std::vector<TraceResult>& results,
                       const std::string& fileName) {
    std::ofstream svg(fileName);
    if (!svg) {
        throw std::runtime_error("Failed to open SVG output file: " + fileName);
    }

    const auto palette = makePalette();
    SvgWriter writer(params,
                     visualization.nearFieldHeight,
                     visualization.nearFieldXMin,
                     visualization.nearFieldXMax,
                     visualization.nearFieldYMin,
                     visualization.nearFieldYMax);
    writer.writeHeader(svg, "Near Field Trace: droplet/background detail");
    writer.writeScene(svg, true);
    for (std::size_t i = 0; i < results.size(); ++i) {
        const std::string& color = palette[i % palette.size()];
        writer.writePath(svg, results[i].noDropletPath, color, true, 1.2);
        writer.writePath(svg, results[i].dropletPath, color, false, 2.0);
    }
    for (std::size_t i = 0; i < results.size(); ++i) {
        writeRefractionAnnotations(svg, writer, results[i], palette[i % palette.size()]);
    }
    writer.writeFooter(svg);
}

void writeCsv(const std::vector<TraceResult>& results, const std::string& fileName) {
    std::ofstream csv(fileName);
    if (!csv) {
        throw std::runtime_error("Failed to open CSV output file: " + fileName);
    }

    csv << "sample_index,sensor_y_mm,status,"
        << "Q0_x_mm,Q0_y_mm,Q_x_mm,Q_y_mm,delta_x_mm,delta_y_mm,"
        << "P1_x_mm,P1_y_mm,P2_x_mm,P2_y_mm\n";
    for (const auto& r : results) {
        csv << r.index << ","
            << formatNumber(r.sensorPoint.y) << ","
            << r.status << ","
            << formatNumber(r.q0.x) << ","
            << formatNumber(r.q0.y) << ","
            << formatNumber(r.q.x) << ","
            << formatNumber(r.q.y) << ","
            << formatNumber(r.delta.x) << ","
            << formatNumber(r.delta.y) << ","
            << formatNumber(r.p1.x) << ","
            << formatNumber(r.p1.y) << ","
            << formatNumber(r.p2.x) << ","
            << formatNumber(r.p2.y) << "\n";
    }
}

TraceSummary summarizeResults(const std::vector<TraceResult>& results) {
    TraceSummary summary;
    summary.total = static_cast<int>(results.size());

    for (const auto& r : results) {
        if (r.status == "no_hit") {
            ++summary.noHit;
        } else if (r.status == "refracted") {
            ++summary.refracted;
        } else if (r.status == "TIR") {
            ++summary.tir;
        } else if (r.status == "invalid") {
            ++summary.invalid;
        }

        const double deltaMagnitude = norm(r.delta);
        summary.minDeltaY = std::min(summary.minDeltaY, r.delta.y);
        summary.maxDeltaY = std::max(summary.maxDeltaY, r.delta.y);
        summary.minDeltaMagnitude = std::min(summary.minDeltaMagnitude, deltaMagnitude);
        summary.maxDeltaMagnitude = std::max(summary.maxDeltaMagnitude, deltaMagnitude);
    }

    if (results.empty()) {
        summary.minDeltaY = 0.0;
        summary.maxDeltaY = 0.0;
        summary.minDeltaMagnitude = 0.0;
        summary.maxDeltaMagnitude = 0.0;
    }

    return summary;
}

void printParameters(const SystemParameters& params) {
    std::cout << "2D coaxial BOS ray tracing parameters (unit: mm)\n";
    std::cout << "  Background plane: x=" << params.backgroundX << "\n";
    std::cout << "  Camera center: (" << params.cameraCenter.x << ", "
              << params.cameraCenter.y << ")\n";
    std::cout << "  Sensor center: (" << params.sensorCenter.x << ", "
              << params.sensorCenter.y << ")\n";
    std::cout << "  Sensor height: " << params.sensorSizeY << "\n";
    std::cout << "  Pixel size: " << params.pixelSizeMm
              << ", sampling step: " << params.sampleStepPixels << " px = "
              << params.pixelSizeMm * params.sampleStepPixels << " mm\n";
    std::cout << "  Droplet center: (" << params.dropletCenter.x << ", "
              << params.dropletCenter.y << "), radius=" << params.dropletRadius << "\n";
    std::cout << "  Refractive index: air=" << params.refractiveIndexAir
              << ", water=" << params.refractiveIndexWater << "\n\n";
}

void printResults(const std::vector<TraceResult>& results) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Trace results\n";
    std::cout << "idx  sensor_y   status       Q0(x,y)                 Q(x,y)"
              << "                  Delta(x,y)\n";
    std::cout << "-------------------------------------------------------------------------------------\n";
    for (const auto& r : results) {
        std::cout << std::setw(3) << r.index << "  "
                  << std::setw(8) << r.sensorPoint.y << "  "
                  << std::left << std::setw(11) << r.status << std::right
                  << "(" << std::setw(10) << r.q0.x << ", " << std::setw(10) << r.q0.y << ")  "
                  << "(" << std::setw(10) << r.q.x << ", " << std::setw(10) << r.q.y << ")  "
                  << "(" << std::setw(10) << r.delta.x << ", " << std::setw(10) << r.delta.y << ")\n";
    }
}

void printSummary(const TraceSummary& summary) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nTrace summary\n";
    std::cout << "  Total samples: " << summary.total << "\n";
    std::cout << "  Status counts: no_hit=" << summary.noHit
              << ", refracted=" << summary.refracted
              << ", TIR=" << summary.tir
              << ", invalid=" << summary.invalid << "\n";
    std::cout << "  DeltaQ.y range: [" << summary.minDeltaY
              << ", " << summary.maxDeltaY << "] mm\n";
    std::cout << "  |DeltaQ| range: [" << summary.minDeltaMagnitude
              << ", " << summary.maxDeltaMagnitude << "] mm\n";
}

void printRefractionDetails(const std::vector<TraceResult>& results) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nRefraction detail\n";
    std::cout << "idx  sensor_y   P1(x,y)                  entry angles(deg)"
              << "        P2(x,y)                  exit angles(deg)\n";
    std::cout << "------------------------------------------------------------------------------------------------\n";

    for (const auto& r : results) {
        if (r.status != "refracted") {
            continue;
        }
        std::cout << std::setw(3) << r.index << "  "
                  << std::setw(8) << r.sensorPoint.y << "  "
                  << "(" << std::setw(10) << r.p1.x << ", " << std::setw(10) << r.p1.y << ")  "
                  << "incident=" << std::setw(9) << r.entryIncidentAngleDeg
                  << ", refracted=" << std::setw(9) << r.entryRefractedAngleDeg << "  "
                  << "(" << std::setw(10) << r.p2.x << ", " << std::setw(10) << r.p2.y << ")  "
                  << "incident=" << std::setw(9) << r.exitIncidentAngleDeg
                  << ", refracted=" << std::setw(9) << r.exitRefractedAngleDeg << "\n";
    }
}

int main() {
    try {
        const SystemParameters params;
        const OutputFiles outputFiles;
        const VisualizationSettings visualization;
        const auto sensorPoints = makeSensorPoints(params);

        std::vector<TraceResult> results;
        results.reserve(sensorPoints.size());
        for (std::size_t i = 0; i < sensorPoints.size(); ++i) {
            results.push_back(traceSample(params, static_cast<int>(i), sensorPoints[i]));
        }

        printParameters(params);
        printResults(results);
        printSummary(summarizeResults(results));
        printRefractionDetails(results);

        writeTraceSvg(params,
                      visualization,
                      results,
                      outputFiles.noDropletSvg,
                      "No Droplet Ray Trace",
                      false,
                      false,
                      false);
        writeTraceSvg(params,
                      visualization,
                      results,
                      outputFiles.dropletSvg,
                      "With Droplet Ray Trace",
                      true,
                      true,
                      false);
        writeComparisonSvg(params, visualization, results, outputFiles.comparisonSvg);
        writeNearFieldSvg(params, visualization, results, outputFiles.nearFieldSvg);
        writeCsv(results, outputFiles.displacementCsv);

        std::cout << "\nWrote outputs:\n";
        std::cout << "  " << outputFiles.noDropletSvg << "\n";
        std::cout << "  " << outputFiles.dropletSvg << "\n";
        std::cout << "  " << outputFiles.comparisonSvg << "\n";
        std::cout << "  " << outputFiles.nearFieldSvg << "\n";
        std::cout << "  " << outputFiles.displacementCsv << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
