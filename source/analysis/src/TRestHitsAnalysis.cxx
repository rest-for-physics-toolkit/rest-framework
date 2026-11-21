#include "TRestHitsAnalysis.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

using TRestHitsAnalysis::XYZVector;

namespace {
using Cell = std::array<long long, 3>;

struct CellHash {
    size_t operator()(const Cell& c) const {
        return std::hash<long long>()(c[0] * 73856093LL ^ c[1] * 19349663LL ^ c[2] * 83492791LL);
    }
};

double GetCoordinate(const XYZVector& p, int axis) { return axis == 0 ? p.X() : axis == 1 ? p.Y() : p.Z(); }

// Cell of side cellSize of a point along the enabled coordinates (0 along the others)
Cell GetCell(const XYZVector& p, double cellSize, const std::array<bool, 3>& coordinates) {
    Cell cell = {0, 0, 0};
    for (int axis = 0; axis < 3; axis++) {
        if (coordinates[axis]) cell[axis] = (long long)std::floor(GetCoordinate(p, axis) / cellSize);
    }
    return cell;
}

// Offsets of a cell and the cells touching it along the enabled coordinates
std::vector<Cell> GetNeighbourOffsets(const std::array<bool, 3>& coordinates) {
    std::vector<Cell> offsets = {{0, 0, 0}};
    for (int axis = 0; axis < 3; axis++) {
        if (!coordinates[axis]) continue;
        std::vector<Cell> expanded;
        for (const auto& offset : offsets) {
            for (int d = -1; d <= 1; d++) {
                Cell neighbour = offset;
                neighbour[axis] = d;
                expanded.push_back(neighbour);
            }
        }
        offsets = std::move(expanded);
    }
    return offsets;
}

double GetPathLength(const std::vector<size_t>& path, const std::vector<std::vector<double>>& distance,
                     bool cyclic) {
    double length = 0;
    for (size_t i = 1; i < path.size(); i++) length += distance[path[i - 1]][path[i]];
    if (cyclic && path.size() > 1) length += distance[path.back()][path.front()];
    return length;
}

// Held-Karp: shortest path visiting the subset `mask` and ending at `j`, built subset by subset. Open
// paths start anywhere, cyclic ones at point 0 (and come back to it)
std::vector<size_t> GetShortestPathExact(const std::vector<std::vector<double>>& distance, bool cyclic) {
    const size_t n = distance.size();
    const size_t full = (size_t)1 << n;
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> cost(full * n, inf);
    std::vector<int> previous(full * n, -1);
    if (cyclic) {
        cost[1 * n + 0] = 0;
    } else {
        for (size_t j = 0; j < n; j++) cost[((size_t)1 << j) * n + j] = 0;
    }

    for (size_t mask = 1; mask < full; mask++) {
        for (size_t j = 0; j < n; j++) {
            const double c = cost[mask * n + j];
            if (c == inf) continue;
            for (size_t k = 0; k < n; k++) {
                if (mask & ((size_t)1 << k)) continue;
                const size_t next = (mask | ((size_t)1 << k)) * n + k;
                if (c + distance[j][k] < cost[next]) {
                    cost[next] = c + distance[j][k];
                    previous[next] = (int)j;
                }
            }
        }
    }

    size_t last = 0;
    double best = inf;
    for (size_t j = 0; j < n; j++) {
        const double c = cost[(full - 1) * n + j] + (cyclic ? distance[j][0] : 0);
        if (c < best) {
            best = c;
            last = j;
        }
    }

    std::vector<size_t> path;
    size_t mask = full - 1;
    int j = (int)last;
    while (j >= 0) {
        path.push_back(j);
        const int p = previous[mask * n + j];
        mask &= ~((size_t)1 << j);
        j = p;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

// Nearest neighbour path from `start`
std::vector<size_t> GetNearestNeighbourPath(const std::vector<std::vector<double>>& distance, size_t start) {
    const size_t n = distance.size();
    std::vector<bool> visited(n, false);
    std::vector<size_t> path = {start};
    visited[start] = true;
    while (path.size() < n) {
        const size_t current = path.back();
        size_t closest = n;
        for (size_t k = 0; k < n; k++) {
            if (!visited[k] && (closest == n || distance[current][k] < distance[current][closest]))
                closest = k;
        }
        visited[closest] = true;
        path.push_back(closest);
    }
    return path;
}

// 2-opt: reverses the segments [i, j] of the path while that makes it shorter. A cyclic path keeps its
// first point
void ImprovePath2Opt(std::vector<size_t>& path, const std::vector<std::vector<double>>& distance,
                     bool cyclic) {
    const size_t n = path.size();
    bool improved = true;
    while (improved) {
        improved = false;
        for (size_t i = cyclic ? 1 : 0; i + 1 < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                // Edges before i and after j, if they exist (the cyclic path closes at the first point)
                const bool hasBefore = i > 0;
                const bool hasAfter = j + 1 < n || cyclic;
                const size_t before = hasBefore ? path[i - 1] : 0;
                const size_t after = j + 1 < n ? path[j + 1] : path[0];
                double delta = 0;
                if (hasBefore) delta += distance[before][path[j]] - distance[before][path[i]];
                if (hasAfter) delta += distance[path[i]][after] - distance[path[j]][after];
                if (delta < -1e-12) {
                    std::reverse(path.begin() + i, path.begin() + j + 1);
                    improved = true;
                }
            }
        }
    }
}

// Line that best fits the points: energy weighted centroid and principal axis of the energy weighted
// covariance (unit vector, null for a single point). False if there are no points or no energy
bool FitLine(const std::vector<XYZVector>& points, const std::vector<double>& energies, XYZVector& centroid,
             XYZVector& axis) {
    if (points.empty() || points.size() != energies.size()) return false;

    double totalEnergy = 0;
    centroid = XYZVector();
    for (size_t i = 0; i < points.size(); i++) {
        const double w = std::max(energies[i], 0.);
        totalEnergy += w;
        centroid += w * points[i];
    }
    if (totalEnergy <= 0) return false;
    centroid /= totalEnergy;

    double covariance[3][3] = {};
    for (size_t i = 0; i < points.size(); i++) {
        const double w = std::max(energies[i], 0.) / totalEnergy;
        const XYZVector d = points[i] - centroid;
        const double c[3] = {d.X(), d.Y(), d.Z()};
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) covariance[a][b] += w * c[a] * c[b];
        }
    }
    // Direction of largest spread: eigenvector of the largest eigenvalue (TMatrixDSymEigen sorts them in
    // descending order), with the sign of the column of the covariance with the largest norm. A single
    // point (no spread) gives no direction
    TMatrixDSym matrix(3);
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) matrix(a, b) = covariance[a][b];
    }
    const TMatrixDSymEigen eigen(matrix);
    const TVectorD& eigenValues = eigen.GetEigenValues();
    const TMatrixD& eigenVectors = eigen.GetEigenVectors();
    axis = XYZVector();
    if (eigenValues[0] > 0) {
        axis = XYZVector(eigenVectors(0, 0), eigenVectors(1, 0), eigenVectors(2, 0)).Unit();
        XYZVector largestColumn;
        for (int a = 0; a < 3; a++) {
            const XYZVector column(covariance[0][a], covariance[1][a], covariance[2][a]);
            if (column.R() > largestColumn.R()) largestColumn = column;
        }
        if (axis.Dot(largestColumn) < 0) axis = -axis;
    }
    return true;
}
}  // namespace

std::vector<std::vector<size_t>> TRestHitsAnalysis::FindClusters(const std::vector<XYZVector>& points,
                                                                 double distance,
                                                                 const std::array<bool, 3>& coordinates) {
    // distance <= 0 (or NaN): no points are joined, every point is a cluster
    if (!(distance > 0)) {
        std::vector<std::vector<size_t>> clusters(points.size());
        for (size_t n = 0; n < points.size(); n++) clusters[n] = {n};
        return clusters;
    }
    const double maxDistance2 = distance * distance;
    auto distance2 = [&](size_t a, size_t b) {
        double d2 = 0;
        for (int axis = 0; axis < 3; axis++) {
            if (!coordinates[axis]) continue;
            const double d = GetCoordinate(points[a], axis) - GetCoordinate(points[b], axis);
            d2 += d * d;
        }
        return d2;
    };

    // Cells of side distance along the enabled coordinates (0 along the others): the points closer than
    // distance to a point are in its cell or in the adjacent ones
    const double cellSize = distance;
    auto getCell = [&](size_t n) { return GetCell(points[n], cellSize, coordinates); };

    std::unordered_map<Cell, std::vector<size_t>, CellHash> cells;
    for (size_t n = 0; n < points.size(); n++) cells[getCell(n)].push_back(n);

    const std::vector<Cell> offsets = GetNeighbourOffsets(coordinates);

    std::vector<std::vector<size_t>> clusters;
    std::vector<bool> assigned(points.size(), false);
    for (size_t seed = 0; seed < points.size(); seed++) {
        if (assigned[seed]) continue;

        // Breadth-first search of the points closer than distance to any point already in the cluster
        std::vector<size_t> cluster = {seed};
        assigned[seed] = true;
        for (size_t q = 0; q < cluster.size(); q++) {
            const size_t point = cluster[q];
            const Cell cell = getCell(point);
            for (const auto& offset : offsets) {
                auto it = cells.find({cell[0] + offset[0], cell[1] + offset[1], cell[2] + offset[2]});
                if (it == cells.end()) continue;

                // Assigned points leave the cell, so they are not compared again
                auto& members = it->second;
                for (size_t m = 0; m < members.size();) {
                    const size_t j = members[m];
                    if (!assigned[j]) {
                        if (distance2(point, j) >= maxDistance2) {
                            m++;
                            continue;
                        }
                        assigned[j] = true;
                        cluster.push_back(j);
                    }
                    members[m] = members.back();
                    members.pop_back();
                }
            }
        }
        std::sort(cluster.begin(), cluster.end());
        clusters.push_back(std::move(cluster));
    }
    return clusters;
}

std::vector<size_t> TRestHitsAnalysis::GetShortestPath(const std::vector<XYZVector>& points, bool cyclic,
                                                       size_t maxExactPoints) {
    const size_t n = points.size();
    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    if (n < 3) return order;

    std::vector<std::vector<double>> distance(n, std::vector<double>(n, 0));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) distance[i][j] = distance[j][i] = (points[i] - points[j]).R();
    }

    if (n <= std::min<size_t>(maxExactPoints, 16)) return GetShortestPathExact(distance, cyclic);

    // Nearest neighbour from every point, or from the farthest point from the centroid
    std::vector<size_t> starts;
    if (n <= 200) {
        starts = order;
    } else {
        XYZVector centroid;
        for (const auto& p : points) centroid += p;
        centroid /= (double)n;
        size_t farthest = 0;
        for (size_t i = 1; i < n; i++) {
            if ((points[i] - centroid).Mag2() > (points[farthest] - centroid).Mag2()) farthest = i;
        }
        starts = {farthest};
    }

    std::vector<size_t> best;
    double bestLength = std::numeric_limits<double>::infinity();
    for (size_t start : starts) {
        std::vector<size_t> path = GetNearestNeighbourPath(distance, start);
        const double length = GetPathLength(path, distance, cyclic);
        if (length < bestLength) {
            bestLength = length;
            best = std::move(path);
        }
    }
    ImprovePath2Opt(best, distance, cyclic);
    return best;
}

std::vector<TRestHitsAnalysis::Node> TRestHitsAnalysis::KMeansClustering(
    const std::vector<XYZVector>& points, const std::vector<double>& energies,
    const std::vector<XYZVector>& centroids, int maxIt, bool fixBoundaries) {
    const int nPoints = (int)points.size();
    const int nodes = (int)centroids.size();
    if (nodes < 2 || nPoints == 0) return {};

    auto isFixed = [&](int n) { return fixBoundaries && (n == 0 || n == nodes - 1); };

    std::vector<XYZVector> centroid = centroids, centroidOld(nodes);

    std::vector<std::vector<int>> clusterPoints(nodes);
    for (int it = 0; it < maxIt; it++) {
        for (auto& cluster : clusterPoints) cluster.clear();

        for (int i = 0; i < nPoints; i++) {
            double minDistance = std::numeric_limits<double>::max();
            int closest = -1;
            for (int n = 0; n < nodes; n++) {
                if (isFixed(n)) continue;
                const double d = (centroid[n] - points[i]).R();
                if (d < minDistance) {
                    minDistance = d;
                    closest = n;
                }
            }
            if (closest >= 0) clusterPoints[closest].push_back(i);
        }

        bool converge = true;
        for (int n = 0; n < nodes; n++) {
            if (isFixed(n)) continue;
            XYZVector mean;
            double totalEnergy = 0;
            for (int i : clusterPoints[n]) {
                mean += energies[i] * points[i];
                totalEnergy += energies[i];
            }
            if (totalEnergy > 0) centroid[n] = mean / totalEnergy;
            converge &= (centroid[n] - centroidOld[n]).R() < 1e-9;
            centroidOld[n] = centroid[n];
        }
        if (converge) break;
    }

    std::vector<Node> result;
    for (int n = 0; n < nodes; n++) {
        if (isFixed(n)) {
            result.push_back({centroid[n], 0, XYZVector()});
        } else if (!clusterPoints[n].empty()) {
            double energy = 0;
            double variance[3] = {};
            for (int i : clusterPoints[n]) {
                energy += energies[i];
                const XYZVector d = points[i] - centroid[n];
                variance[0] += energies[i] * d.X() * d.X();
                variance[1] += energies[i] * d.Y() * d.Y();
                variance[2] += energies[i] * d.Z() * d.Z();
            }
            auto sigma = [&](int a) {
                return energy > 0 ? std::sqrt(std::max(variance[a] / energy, 0.)) : 0;
            };
            result.push_back({centroid[n], energy, XYZVector(sigma(0), sigma(1), sigma(2))});
        }
    }
    return result;
}

std::vector<TRestHitsAnalysis::Node> TRestHitsAnalysis::MergePoints(
    const std::vector<XYZVector>& points, const std::vector<double>& energies, double startingDistance,
    double minimumDistance, double distanceStepFactor, size_t maxNodes,
    const std::vector<XYZVector>& sigmas) {
    if (startingDistance <= 0 || distanceStepFactor <= 1) {
        throw std::invalid_argument(
            "TRestHitsAnalysis::MergePoints: startingDistance must be > 0 and distanceStepFactor > 1");
    }
    if (points.empty() || points.size() != energies.size()) return {};
    maxNodes = std::max<size_t>(maxNodes, 1);

    struct Group {
        XYZVector position;
        double energy;
        std::vector<size_t> members;
    };
    std::vector<Group> groups;
    groups.reserve(points.size());
    for (size_t i = 0; i < points.size(); i++) groups.push_back({points[i], energies[i], {i}});

    // As legacy TRestHits::MergeHits: energy weighted position, the first one kept without energy
    auto merge = [](Group& a, Group& b) {
        const double totalEnergy = a.energy + b.energy;
        if (totalEnergy > 0) a.position = (a.energy * a.position + b.energy * b.position) / totalEnergy;
        a.energy = totalEnergy;
        a.members.insert(a.members.end(), b.members.begin(), b.members.end());
    };

    double distance = startingDistance;
    while (distance < minimumDistance || groups.size() > maxNodes) {
        const double distance2 = distance * distance;
        bool merged = true;
        while (merged) {
            merged = false;
            for (size_t i = 0; i < groups.size(); i++) {
                for (size_t j = i + 1; j < groups.size();) {
                    if ((groups[i].position - groups[j].position).Mag2() < distance2) {
                        merge(groups[i], groups[j]);
                        groups.erase(groups.begin() + j);
                        merged = true;
                    } else {
                        j++;
                    }
                }
            }
        }
        distance *= distanceStepFactor;
    }

    std::vector<Node> nodes;
    nodes.reserve(groups.size());
    for (const auto& group : groups) {
        double variance[3] = {};
        for (size_t i : group.members) {
            const XYZVector d = points[i] - group.position;
            const XYZVector s = i < sigmas.size() ? sigmas[i] : XYZVector();
            variance[0] += energies[i] * (d.X() * d.X() + s.X() * s.X());
            variance[1] += energies[i] * (d.Y() * d.Y() + s.Y() * s.Y());
            variance[2] += energies[i] * (d.Z() * d.Z() + s.Z() * s.Z());
        }
        auto sigma = [&](int a) {
            return group.energy > 0 ? std::sqrt(std::max(variance[a] / group.energy, 0.)) : 0;
        };
        nodes.push_back({group.position, group.energy, XYZVector(sigma(0), sigma(1), sigma(2))});
    }
    return nodes;
}

bool TRestHitsAnalysis::GetLine(const std::vector<XYZVector>& points, const std::vector<double>& energies,
                                XYZVector& first, XYZVector& last) {
    XYZVector centroid, axis;
    if (!FitLine(points, energies, centroid, axis)) return false;

    double tMin = 0, tMax = 0;
    for (size_t i = 0; i < points.size(); i++) {
        const double t = (points[i] - centroid).Dot(axis);
        if (i == 0 || t < tMin) tMin = t;
        if (i == 0 || t > tMax) tMax = t;
    }
    first = centroid + tMin * axis;
    last = centroid + tMax * axis;
    return true;
}

std::vector<XYZVector> TRestHitsAnalysis::GetLineNodes(const std::vector<XYZVector>& points,
                                                       const std::vector<double>& energies, int nodes) {
    XYZVector first, last;
    if (nodes < 2 || !GetLine(points, energies, first, last)) return {};

    std::vector<XYZVector> result(nodes);
    for (int n = 0; n < nodes; n++) result[n] = first + (double)n / (nodes - 1) * (last - first);
    return result;
}

double TRestHitsAnalysis::GetLineResidual(const std::vector<XYZVector>& points,
                                          const std::vector<double>& energies) {
    XYZVector centroid, axis;
    if (!FitLine(points, energies, centroid, axis)) return -1;

    double sum = 0, totalEnergy = 0;
    for (size_t i = 0; i < points.size(); i++) {
        const double w = std::max(energies[i], 0.);
        const XYZVector d = points[i] - centroid;
        const XYZVector perpendicular = d - d.Dot(axis) * axis;
        sum += w * perpendicular.Mag2();
        totalEnergy += w;
    }
    return std::sqrt(sum / totalEnergy);
}

double TRestHitsAnalysis::GetTrackEnds(const std::vector<XYZVector>& points,
                                       const std::vector<double>& energies, XYZVector& origin,
                                       XYZVector& end) {
    if (points.empty()) return 0;
    const XYZVector& first = points.front();
    const XYZVector& last = points.back();
    const XYZVector segment = last - first;
    const double length2 = segment.Mag2();

    // Energy of the halves closer to the first and to the last point
    double firstEnergy = 0, lastEnergy = 0;
    for (size_t i = 0; i < points.size(); i++) {
        const double fraction = length2 > 0 ? (points[i] - first).Dot(segment) / length2
                                            : (points.size() > 1 ? (double)i / (points.size() - 1) : 0);
        if (fraction < 0.5) {
            firstEnergy += energies[i];
        } else if (fraction > 0.5) {
            lastEnergy += energies[i];
        } else {
            firstEnergy += energies[i] / 2.;
            lastEnergy += energies[i] / 2.;
        }
    }

    if (firstEnergy > lastEnergy) {
        end = first;
        origin = last;
    } else {
        origin = first;
        end = last;
    }
    const double total = firstEnergy + lastEnergy;
    return total > 0 ? std::abs(lastEnergy - firstEnergy) / total : 0;
}

namespace {
/// Minimum of f in [low, high], deterministic: the best of nScan + 1 equidistant points, then a golden
/// section search between its neighbours. Returns {x, f(x)} and the index of the best point of the scan
std::pair<double, double> MinimizeOnInterval(const std::function<double(double)>& f, double low, double high,
                                             int nScan, int nGolden, int& best) {
    const double step = (high - low) / nScan;
    best = 0;
    double bestValue = std::numeric_limits<double>::max();
    for (int i = 0; i <= nScan; i++) {
        const double value = f(low + i * step);
        if (value < bestValue) {
            bestValue = value;
            best = i;
        }
    }
    const double ratio = (std::sqrt(5.0) - 1) / 2;
    double a = low + std::max(best - 1, 0) * step;
    double b = low + std::min(best + 1, nScan) * step;
    double c = b - ratio * (b - a);
    double d = a + ratio * (b - a);
    double fC = f(c);
    double fD = f(d);
    for (int i = 0; i < nGolden; i++) {
        if (fC < fD) {
            b = d;
            d = c;
            fD = fC;
            c = b - ratio * (b - a);
            fC = f(c);
        } else {
            a = c;
            c = d;
            fC = fD;
            d = a + ratio * (b - a);
            fD = f(d);
        }
    }
    const double x = (a + b) / 2;
    return {x, f(x)};
}
}  // namespace

double TRestHitsAnalysis::GetGaussSigma(const std::vector<double>& values,
                                        const std::vector<double>& energies) {
    if (values.empty()) return 0;

    // Energy of each position (several hits in a strip are added); positions closer than a millionth
    // of their size are the same one
    std::map<double, double> byPosition;
    for (size_t n = 0; n < values.size() && n < energies.size(); n++) byPosition[values[n]] += energies[n];
    std::vector<std::pair<double, double>> positions;
    for (const auto& [x, e] : byPosition) {
        if (!positions.empty() && x - positions.back().first <= 1e-6 * std::max(1.0, std::abs(x))) {
            positions.back().second += e;
        } else {
            positions.emplace_back(x, e);
        }
    }

    double total = 0;
    double pitch = std::numeric_limits<double>::max();
    for (size_t n = 0; n < positions.size(); n++) {
        total += positions[n].second;
        if (n > 0) pitch = std::min(pitch, positions[n].first - positions[n - 1].first);
    }
    if (!(total > 0)) return -1;
    // With two strips any sigma below the pitch fits (the centre moving with it): three are needed
    const auto significant = std::count_if(positions.begin(), positions.end(),
                                           [total](const auto& p) { return p.second > 0.01 * total; });
    if (significant < 3) return -1;
    const double xMin = positions.front().first;
    const double xMax = positions.back().first;

    // Strips of width pitch from one before the first to one after the last (0 where there is no hit)
    const long nStrips = std::lround((xMax - xMin) / pitch) + 3;
    if (nStrips > 100000) return -1;
    std::vector<double> stripEnergy(nStrips, 0);
    for (const auto& [x, e] : positions) stripEnergy[std::lround((x - xMin) / pitch) + 1] += e;

    // Least squares of the energies against the gaussian of integral total integrated over each strip
    const double halfPitch = pitch / 2;
    const auto chi2 = [&](double mean, double sigma) {
        const double scale = 1 / (sigma * std::sqrt(2.0));
        double sum = 0;
        for (long k = 0; k < nStrips; k++) {
            const double x = xMin + (k - 1) * pitch - mean;
            const double expected =
                total * 0.5 * (std::erf((x + halfPitch) * scale) - std::erf((x - halfPitch) * scale));
            sum += (stripEnergy[k] - expected) * (stripEnergy[k] - expected);
        }
        return sum;
    };

    // For each sigma (in log, from pitch / 100 to the extent) the best centre, within the strips
    const double extent = xMax - xMin + 2 * pitch;
    const auto bestForSigma = [&](double logSigma) {
        const double sigma = std::exp(logSigma);
        int best = 0;
        return MinimizeOnInterval([&](double mean) { return chi2(mean, sigma); }, xMin - halfPitch,
                                  xMax + halfPitch, 40, 40, best)
            .second;
    };
    constexpr int nScan = 40;
    int best = 0;
    const double logSigma =
        MinimizeOnInterval(bestForSigma, std::log(pitch / 100), std::log(extent), nScan, 40, best).first;
    // At the extent: the energy does not decrease away from the centre (no gaussian)
    if (best == nScan) return -1;
    return std::exp(logSigma);
}

double TRestHitsAnalysis::GetTwist(const std::vector<XYZVector>& points, size_t n, size_t m) {
    if (m == 0) m = points.size();
    if (m < n + 2) return 0;

    double sum = 0;
    int count = 0;
    for (size_t N = n + 1; N + 1 < m; N++) {
        // Unit: null vector if the points coincide (as TVector3::Unit in legacy)
        const XYZVector a = (points[N - 1] - points[N + 1]).Unit();
        const XYZVector b = (points[N] - points[N + 1]).Unit();
        sum += (1 - a.Dot(b)) / 2.;
        count++;
    }
    return count == 0 ? -1 : sum / count;
}

double TRestHitsAnalysis::GetTwistWeighted(const std::vector<XYZVector>& points,
                                           const std::vector<double>& energies, size_t n, size_t m) {
    if (m == 0) m = points.size();
    if (m < n + 2) return 0;

    double sum = 0;
    int count = 0;
    for (size_t N = n + 1; N + 1 < m; N++) {
        const XYZVector a = (points[N - 1] - points[N]).Unit();
        const XYZVector b = (points[N] - points[N + 1]).Unit();
        const double weight =
            std::abs(energies[N - 1] - energies[N]) + std::abs(energies[N] - energies[N + 1]);
        sum += weight * (1 - a.Dot(b)) / 2.;
        count++;
    }
    return count == 0 ? -1 : sum / count;
}
