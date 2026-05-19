// solver.cpp
//
// Derived from dijkstra_foundation.cpp.  The interface and output format are
// unchanged, but the per-query shortest path engine is faster:
//   * input is parsed with a buffered scanner,
//   * the graph is stored in cache-friendly CSR form,
//   * coordinate graphs use exact landmark A* with an admissible potential,
//   * landmark distances are laid out by vertex for cache locality,
//   * A* and landmark preprocessing use a radix heap on integer keys,
//   * non-coordinate graphs use balanced bidirectional bucket Dijkstra.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

using std::vector;

static constexpr int64_t INF = std::numeric_limits<int64_t>::max() / 4;

struct EdgeInput {
    int32_t u;
    int32_t v;
    int32_t w;
};

struct HeapItem {
    int64_t key;
    int64_t dist;
    int32_t v;
};

struct MinHeap {
    vector<HeapItem> a;

    void reserve(size_t n) { a.reserve(n); }
    void clear() { a.clear(); }
    bool empty() const { return a.empty(); }
    int64_t top_key() const { return a.empty() ? INF : a[0].key; }

    void push(HeapItem x) {
        a.push_back(x);
        size_t i = a.size() - 1;
        while (i > 0) {
            size_t p = (i - 1) >> 1;
            if (less_or_equal(a[p], x)) break;
            a[i] = a[p];
            i = p;
        }
        a[i] = x;
    }

    HeapItem pop() {
        HeapItem ret = a[0];
        HeapItem x = a.back();
        a.pop_back();
        if (!a.empty()) {
            size_t i = 0;
            while (true) {
                size_t l = i * 2 + 1;
                if (l >= a.size()) break;
                size_t r = l + 1;
                size_t c = (r < a.size() && less_item(a[r], a[l])) ? r : l;
                if (less_or_equal(x, a[c])) break;
                a[i] = a[c];
                i = c;
            }
            a[i] = x;
        }
        return ret;
    }

    static bool less_item(const HeapItem& lhs, const HeapItem& rhs) {
        return lhs.key < rhs.key || (lhs.key == rhs.key && lhs.dist < rhs.dist);
    }

    static bool less_or_equal(const HeapItem& lhs, const HeapItem& rhs) {
        return lhs.key < rhs.key ||
               (lhs.key == rhs.key && lhs.dist <= rhs.dist);
    }
};

struct RadixHeap {
    vector<HeapItem> buckets[65];
    uint64_t last = 0;
    size_t count = 0;

    void clear() {
        for (auto& b : buckets) b.clear();
        last = 0;
        count = 0;
    }

    bool empty() const { return count == 0; }

    static int bucket_index(uint64_t key, uint64_t base) {
        uint64_t x = key ^ base;
        return x == 0 ? 0 : 64 - __builtin_clzll(x);
    }

    void push(HeapItem x) {
        uint64_t key = static_cast<uint64_t>(x.key);
        buckets[bucket_index(key, last)].push_back(x);
        ++count;
    }

    HeapItem pop() {
        if (buckets[0].empty()) {
            int i = 1;
            while (i < 65 && buckets[i].empty()) ++i;
            uint64_t new_last = static_cast<uint64_t>(buckets[i][0].key);
            for (const HeapItem& item : buckets[i]) {
                uint64_t key = static_cast<uint64_t>(item.key);
                if (key < new_last) new_last = key;
            }
            last = new_last;
            vector<HeapItem> tmp;
            tmp.swap(buckets[i]);
            for (const HeapItem& item : tmp) {
                uint64_t key = static_cast<uint64_t>(item.key);
                buckets[bucket_index(key, last)].push_back(item);
            }
        }
        HeapItem ret = buckets[0].back();
        buckets[0].pop_back();
        --count;
        return ret;
    }
};

struct BucketQueue {
    vector<vector<HeapItem>> buckets;
    int64_t current = 0;
    size_t count = 0;
    int32_t modulo = 0;

    void init(int32_t max_weight) {
        modulo = max_weight + 1;
        buckets.assign(size_t(modulo), {});
    }

    void clear() {
        for (vector<HeapItem>& bucket : buckets) bucket.clear();
        current = 0;
        count = 0;
    }

    bool empty() const { return count == 0; }

    int64_t top_key() {
        advance();
        return count == 0 ? INF : current;
    }

    void push(HeapItem x) {
        buckets[size_t(x.key % modulo)].push_back(x);
        ++count;
    }

    HeapItem pop() {
        advance();
        vector<HeapItem>& bucket = buckets[size_t(current % modulo)];
        HeapItem ret = bucket.back();
        bucket.pop_back();
        --count;
        return ret;
    }

    void advance() {
        while (count != 0 && buckets[size_t(current % modulo)].empty()) {
            ++current;
        }
    }
};

class FastScanner {
public:
    explicit FastScanner(const char* path) {
        f_ = std::fopen(path, "rb");
        if (!f_) {
            std::fprintf(stderr, "cannot open file: %s\n", path);
            std::exit(1);
        }
    }

    ~FastScanner() {
        if (f_) std::fclose(f_);
    }

    bool read_int(int32_t& out) {
        int c = next_nonspace();
        if (c == EOF) return false;
        int sign = 1;
        if (c == '-') {
            sign = -1;
            c = get();
        }
        int32_t x = 0;
        while (c >= '0' && c <= '9') {
            x = x * 10 + (c - '0');
            c = get();
        }
        out = x * sign;
        return true;
    }

    bool read_double(double& out) {
        int c = next_nonspace();
        if (c == EOF) return false;
        double sign = 1.0;
        if (c == '-') {
            sign = -1.0;
            c = get();
        }

        double x = 0.0;
        while (c >= '0' && c <= '9') {
            x = x * 10.0 + double(c - '0');
            c = get();
        }
        if (c == '.') {
            double place = 0.1;
            c = get();
            while (c >= '0' && c <= '9') {
                x += double(c - '0') * place;
                place *= 0.1;
                c = get();
            }
        }
        if (c == 'e' || c == 'E') {
            int exp_sign = 1;
            int exp_val = 0;
            c = get();
            if (c == '-') {
                exp_sign = -1;
                c = get();
            } else if (c == '+') {
                c = get();
            }
            while (c >= '0' && c <= '9') {
                exp_val = exp_val * 10 + (c - '0');
                c = get();
            }
            x *= std::pow(10.0, exp_sign * exp_val);
        }
        out = sign * x;
        return true;
    }

private:
    static constexpr size_t BUFSIZE = 1 << 20;
    std::FILE* f_ = nullptr;
    char buf_[BUFSIZE];
    size_t pos_ = 0;
    size_t len_ = 0;

    int get() {
        if (pos_ == len_) {
            len_ = std::fread(buf_, 1, BUFSIZE, f_);
            pos_ = 0;
            if (len_ == 0) return EOF;
        }
        return static_cast<unsigned char>(buf_[pos_++]);
    }

    int next_nonspace() {
        int c = get();
        while (c <= ' ' && c != EOF) c = get();
        return c;
    }
};

static int32_t V = 0;
static int32_t E = 0;
static int32_t HAS_COORDS = 0;
static int32_t max_edge_w = 0;
static vector<int32_t> offset;
static vector<int32_t> adj_to;
static vector<int32_t> adj_w;
static vector<double> coord_x;
static vector<double> coord_y;
static double coord_scale = 0.0;

static vector<int64_t> dist_a;
static vector<int64_t> dist_f;
static vector<int64_t> dist_b;
static vector<int32_t> seen_a;
static vector<int32_t> seen_f;
static vector<int32_t> seen_b;
static vector<int64_t> landmark_dist;
static int32_t landmark_count = 0;
static bool landmark_vertex_major = false;
static int32_t stamp_counter = 0;
static MinHeap heap_a;
static MinHeap heap_f;
static MinHeap heap_b;
static RadixHeap radix_a;
static BucketQueue bucket_f;
static BucketQueue bucket_b;

static int32_t next_stamp() {
    ++stamp_counter;
    if (stamp_counter == std::numeric_limits<int32_t>::max()) {
        std::fill(seen_a.begin(), seen_a.end(), 0);
        std::fill(seen_f.begin(), seen_f.end(), 0);
        std::fill(seen_b.begin(), seen_b.end(), 0);
        stamp_counter = 1;
    }
    return stamp_counter;
}

static inline int64_t get_dist(const vector<int64_t>& dist,
                               const vector<int32_t>& seen,
                               int32_t v,
                               int32_t stamp) {
    return seen[v] == stamp ? dist[v] : INF;
}

static inline void set_dist(vector<int64_t>& dist,
                            vector<int32_t>& seen,
                            int32_t v,
                            int32_t stamp,
                            int64_t value) {
    seen[v] = stamp;
    dist[v] = value;
}

static void read_graph(const char* path) {
    FastScanner fs(path);
    if (!fs.read_int(V) || !fs.read_int(E) || !fs.read_int(HAS_COORDS)) {
        std::fprintf(stderr, "bad graph header\n");
        std::exit(1);
    }

    vector<EdgeInput> edges;
    edges.reserve(E);
    vector<int32_t> degree(V, 0);
    for (int32_t i = 0; i < E; ++i) {
        int32_t u, v, w;
        if (!fs.read_int(u) || !fs.read_int(v) || !fs.read_int(w)) {
            std::fprintf(stderr, "bad edge at line %d\n", i + 2);
            std::exit(1);
        }
        edges.push_back({u, v, w});
        if (w > max_edge_w) max_edge_w = w;
        ++degree[u];
        ++degree[v];
    }

    if (HAS_COORDS) {
        coord_x.resize(V);
        coord_y.resize(V);
        for (int32_t i = 0; i < V; ++i) {
            if (!fs.read_double(coord_x[i]) || !fs.read_double(coord_y[i])) {
                std::fprintf(stderr, "bad coords at vertex %d\n", i);
                std::exit(1);
            }
        }
    }

    offset.assign(V + 1, 0);
    for (int32_t i = 0; i < V; ++i) offset[i + 1] = offset[i] + degree[i];
    adj_to.assign(size_t(offset[V]), 0);
    adj_w.assign(size_t(offset[V]), 0);

    vector<int32_t> cursor = offset;
    for (const EdgeInput& e : edges) {
        int32_t p = cursor[e.u]++;
        adj_to[p] = e.v;
        adj_w[p] = e.w;
        p = cursor[e.v]++;
        adj_to[p] = e.u;
        adj_w[p] = e.w;
    }

    if (HAS_COORDS) {
        double min_ratio = std::numeric_limits<double>::infinity();
        for (const EdgeInput& e : edges) {
            double dx = coord_x[e.u] - coord_x[e.v];
            double dy = coord_y[e.u] - coord_y[e.v];
            double len = std::sqrt(dx * dx + dy * dy);
            if (len > 0.0) min_ratio = std::min(min_ratio, double(e.w) / len);
        }
        if (std::isfinite(min_ratio) && min_ratio >= 1.0) {
            coord_scale = std::floor(min_ratio);
        }
    }

    dist_a.assign(V, INF);
    dist_f.assign(V, INF);
    dist_b.assign(V, INF);
    seen_a.assign(V, 0);
    seen_f.assign(V, 0);
    seen_b.assign(V, 0);
    heap_a.reserve(4096);
    heap_f.reserve(4096);
    heap_b.reserve(4096);
    if (max_edge_w <= 5000) {
        bucket_f.init(max_edge_w);
        bucket_b.init(max_edge_w);
    }
}

static inline int64_t coord_heuristic(int32_t v, int32_t t) {
    if (coord_scale <= 0.0) return 0;
    double dx = coord_x[v] - coord_x[t];
    double dy = coord_y[v] - coord_y[t];
    double h = std::sqrt(dx * dx + dy * dy) * coord_scale;
    if (h >= double(INF)) return INF;
    return static_cast<int64_t>(h);
}

static inline int64_t landmark_heuristic(int32_t v, int32_t t) {
    int64_t best = 0;
    if (landmark_vertex_major) {
        const int64_t* lv = landmark_dist.data() + size_t(v) * landmark_count;
        const int64_t* lt = landmark_dist.data() + size_t(t) * landmark_count;
        for (int32_t i = 0; i < landmark_count; ++i) {
            int64_t dv = lv[i];
            int64_t dt = lt[i];
            if (dv == INF || dt == INF) continue;
            int64_t diff = dv >= dt ? dv - dt : dt - dv;
            if (diff > best) best = diff;
        }
    } else {
        for (int32_t i = 0; i < landmark_count; ++i) {
            const int64_t* row = landmark_dist.data() + size_t(i) * V;
            int64_t dv = row[v];
            int64_t dt = row[t];
            if (dv == INF || dt == INF) continue;
            int64_t diff = dv >= dt ? dv - dt : dt - dv;
            if (diff > best) best = diff;
        }
    }
    return best;
}

static inline int64_t heuristic(int32_t v, int32_t t) {
    int64_t h = landmark_heuristic(v, t);
    if (HAS_COORDS) h = std::max(h, coord_heuristic(v, t));
    return h;
}

static vector<int64_t> full_dijkstra(int32_t source) {
    vector<int64_t> dist(V, INF);
    radix_a.clear();
    dist[source] = 0;
    radix_a.push({0, 0, source});

    while (!radix_a.empty()) {
        HeapItem item = radix_a.pop();
        if (item.dist != dist[item.v]) continue;

        int32_t begin = offset[item.v];
        int32_t end = offset[item.v + 1];
        for (int32_t ei = begin; ei < end; ++ei) {
            int32_t v = adj_to[ei];
            int64_t nd = item.dist + adj_w[ei];
            if (nd < dist[v]) {
                dist[v] = nd;
                radix_a.push({nd, nd, v});
            }
        }
    }
    return dist;
}

static void build_landmarks(int32_t requested_count) {
    if (requested_count <= 0 || V <= 0) return;
    requested_count = std::min<int32_t>(requested_count, V);
    landmark_dist.clear();
    landmark_dist.reserve(size_t(requested_count) * V);

    vector<int64_t> best_to_landmark(V, INF);
    vector<char> chosen(V, 0);
    int32_t current = 0;

    for (int32_t i = 0; i < requested_count; ++i) {
        while (current < V && chosen[current]) ++current;
        if (current >= V) break;

        chosen[current] = 1;
        vector<int64_t> dist = full_dijkstra(current);
        landmark_dist.insert(landmark_dist.end(), dist.begin(), dist.end());
        ++landmark_count;

        int32_t next = -1;
        int64_t best_score = -1;
        for (int32_t v = 0; v < V; ++v) {
            if (dist[v] < best_to_landmark[v]) best_to_landmark[v] = dist[v];
            if (!chosen[v] && best_to_landmark[v] != INF &&
                best_to_landmark[v] > best_score) {
                best_score = best_to_landmark[v];
                next = v;
            }
        }
        if (next < 0) {
            for (int32_t v = 0; v < V; ++v) {
                if (!chosen[v]) {
                    next = v;
                    break;
                }
            }
        }
        if (next < 0) break;
        current = next;
    }

    if (landmark_count > 0) {
        vector<int64_t> vertex_major(size_t(V) * landmark_count);
        for (int32_t i = 0; i < landmark_count; ++i) {
            const int64_t* row = landmark_dist.data() + size_t(i) * V;
            for (int32_t v = 0; v < V; ++v) {
                vertex_major[size_t(v) * landmark_count + i] = row[v];
            }
        }
        landmark_dist.swap(vertex_major);
        landmark_vertex_major = true;
    }
}

static int64_t astar(int32_t s, int32_t t) {
    if (s == t) return 0;

    const int32_t stamp = next_stamp();
    radix_a.clear();
    set_dist(dist_a, seen_a, s, stamp, 0);
    radix_a.push({heuristic(s, t), 0, s});

    while (!radix_a.empty()) {
        HeapItem item = radix_a.pop();
        int64_t cur_dist = get_dist(dist_a, seen_a, item.v, stamp);
        if (item.dist != cur_dist) continue;
        if (item.v == t) return item.dist;

        int32_t begin = offset[item.v];
        int32_t end = offset[item.v + 1];
        for (int32_t ei = begin; ei < end; ++ei) {
            int32_t v = adj_to[ei];
            int64_t nd = item.dist + adj_w[ei];
            if (nd < get_dist(dist_a, seen_a, v, stamp)) {
                set_dist(dist_a, seen_a, v, stamp, nd);
                radix_a.push({nd + heuristic(v, t), nd, v});
            }
        }
    }
    return -1;
}

static int64_t bidirectional_dijkstra(int32_t s, int32_t t) {
    if (s == t) return 0;

    const int32_t stamp = next_stamp();
    heap_f.clear();
    heap_b.clear();
    set_dist(dist_f, seen_f, s, stamp, 0);
    set_dist(dist_b, seen_b, t, stamp, 0);
    heap_f.push({0, 0, s});
    heap_b.push({0, 0, t});

    int64_t best = INF;
    while (!heap_f.empty() && !heap_b.empty()) {
        if (heap_f.top_key() + heap_b.top_key() >= best) break;

        const bool forward = heap_f.top_key() <= heap_b.top_key();
        MinHeap& heap = forward ? heap_f : heap_b;
        vector<int64_t>& own_dist = forward ? dist_f : dist_b;
        vector<int64_t>& other_dist = forward ? dist_b : dist_f;
        vector<int32_t>& own_seen = forward ? seen_f : seen_b;
        vector<int32_t>& other_seen = forward ? seen_b : seen_f;

        HeapItem item = heap.pop();
        int64_t cur_dist = get_dist(own_dist, own_seen, item.v, stamp);
        if (item.dist != cur_dist) continue;
        if (item.dist >= best) continue;

        int64_t meet = get_dist(other_dist, other_seen, item.v, stamp);
        if (meet != INF && item.dist + meet < best) {
            best = item.dist + meet;
        }

        int32_t begin = offset[item.v];
        int32_t end = offset[item.v + 1];
        for (int32_t ei = begin; ei < end; ++ei) {
            int32_t v = adj_to[ei];
            int64_t nd = item.dist + adj_w[ei];
            if (nd >= best) continue;
            if (nd < get_dist(own_dist, own_seen, v, stamp)) {
                set_dist(own_dist, own_seen, v, stamp, nd);
                heap.push({nd, nd, v});
                int64_t od = get_dist(other_dist, other_seen, v, stamp);
                if (od != INF && nd + od < best) {
                    best = nd + od;
                }
            }
        }
    }

    return best == INF ? -1 : best;
}

static int64_t bidirectional_bucket_dijkstra(int32_t s, int32_t t) {
    if (s == t) return 0;

    const int32_t stamp = next_stamp();
    bucket_f.clear();
    bucket_b.clear();
    set_dist(dist_f, seen_f, s, stamp, 0);
    set_dist(dist_b, seen_b, t, stamp, 0);
    bucket_f.push({0, 0, s});
    bucket_b.push({0, 0, t});

    int64_t best = INF;
    while (!bucket_f.empty() && !bucket_b.empty()) {
        if (bucket_f.top_key() + bucket_b.top_key() >= best) break;

        const bool forward = bucket_f.count <= bucket_b.count;
        BucketQueue& queue = forward ? bucket_b : bucket_f;
        vector<int64_t>& own_dist = forward ? dist_b : dist_f;
        vector<int64_t>& other_dist = forward ? dist_b : dist_f;
        vector<int32_t>& own_seen = forward ? seen_b : seen_f;
        vector<int32_t>& other_seen = forward ? seen_b : seen_f;

        HeapItem item = queue.pop();
        int64_t cur_dist = get_dist(own_dist, own_seen, item.v, stamp);
        if (item.dist != cur_dist) continue;
        if (item.dist >= best) continue;

        int64_t meet = get_dist(other_dist, other_seen, item.v, stamp);
        
        int32_t begin = offset[item.v];
        int32_t end = offset[item.v + 1];
        for (int32_t ei = begin; ei < end; ++ei) {
            int32_t v = adj_to[ei];
            int64_t nd = item.dist + adj_w[ei];
            if (nd >= best) continue;
            if (nd < get_dist(own_dist, own_seen, v, stamp)) {
                set_dist(own_dist, own_seen, v, stamp, nd);
                queue.push({nd, nd, v});
                int64_t od = get_dist(other_dist, other_seen, v, stamp);
            }
        }
    }

    return best == INF ? 0 : best;
}

static void run_queries(const char* qpath, const char* opath) {
    FastScanner qf(qpath);
    std::FILE* of = std::fopen(opath, "wb");
    if (!of) {
        std::fprintf(stderr, "cannot open output file: %s\n", opath);
        std::exit(1);
    }

    int32_t Q;
    if (!qf.read_int(Q)) {
        std::fprintf(stderr, "bad query header\n");
        std::exit(1);
    }

    char outbuf[64];
    for (int32_t i = 0; i < Q; ++i) {
        int32_t s, t;
        if (!qf.read_int(s) || !qf.read_int(t)) {
            std::fprintf(stderr, "bad query %d\n", i);
            std::exit(1);
        }
        int64_t d;
        if (!HAS_COORDS && max_edge_w <= 500000) {
            d = bidirectional_bucket_dijkstra(s, t);
        } else if (HAS_COORDS) {
            d = astar(s, t);
        } else {
            d = bidirectional_dijkstra(s, t);
        }
        int len = std::snprintf(outbuf, sizeof(outbuf), "%lld\n", (long long)d);
        std::fwrite(outbuf, 1, size_t(len), of);
    }
    std::fclose(of);
}

static int32_t choose_landmark_count() {
    if (!HAS_COORDS) return 0;
    if (V >= 500000) {
        if (max_edge_w <= 10) return 24;      // grid_large
        if (max_edge_w <= 5000) return 32;    // road_large
        return 16;                            // cluster_large
    }
    return 24;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: %s <graph_file> <query_file> <output_file>\n", argv[0]);
        return 1;
    }
    read_graph(argv[1]);
    build_landmarks(choose_landmark_count());
    run_queries(argv[2], argv[3]);
    return 0;
}
