#include <bits/stdc++.h>

using namespace std;

struct encoder {
    size_t n, k;
    vector<vector<uint64_t>> g;

    encoder(size_t n, size_t k, const vector<vector<uint64_t>> &g_cols)
            : n(n), k(k), g(g_cols) {}

    uint64_t encode(uint64_t word) {
        uint64_t result = 0;
        for (size_t i = 0; i < n; ++i) {
            uint64_t bit = 0;
            for (size_t j = 0; j < k; ++j) {
                bit ^= (word >> j) & g[j][i];
            }
            result |= (bit << i);
        }
        return result;
    }
};

struct trellis_node {
    ssize_t from0 = -1, from1 = -1;
};

struct coder {
    size_t n, k;
    vector<vector<uint64_t>> g;
    vector<vector<uint64_t>> replacement;
    vector<size_t> active_start;
    vector<size_t> active_end;
    vector<vector<trellis_node>> trellis;
    vector<vector<double>> d;
    vector<vector<ssize_t>> prev;

    coder(size_t n, size_t k, const vector<vector<uint64_t>> &g) : n(n), k(k), g(g),
                                                                   replacement(k,vector<uint64_t>(k)),
                                                                   active_start(k), active_end(k),
                                                                   trellis(n + 1), d(n + 1),
                                                                   prev(n + 1) {
        g_msf();
        build_trellis();
    }

    vector<uint64_t> trellis_profile() {
        uint64_t nodes_cnt = 1;
        vector<uint64_t> profile(n + 1);
        profile[0] = 1;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < k; ++j) {
                if (active_start[j] == i) {
                    nodes_cnt <<= 1;
                }
                if (active_end[j] == i) {
                    nodes_cnt >>= 1;
                }
            }
            profile[i + 1] = nodes_cnt;
        }
        return profile;
    }

    uint64_t encode(uint64_t word) {
        uint64_t result = 0;
        for (size_t i = 0; i < n; ++i) {
            uint64_t bit = 0;
            for (size_t j = 0; j < k; ++j) {
                bit ^= (word >> j) & g[j][i];
            }
            result |= (bit << i);
        }
        return result;
    }

    uint64_t decode(vector<double> code, bool flip = false) {
        d[0][0] = 0;
        prev[0][0] = -1;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < trellis[i + 1].size(); ++j) {
                if (trellis[i + 1][j].from0 == -1) {
                    d[i + 1][j] = d[i][trellis[i + 1][j].from1] + code[i];
                    prev[i + 1][j] = trellis[i + 1][j].from1;
                } else if (trellis[i + 1][j].from1 == -1) {
                    d[i + 1][j] = d[i][trellis[i + 1][j].from0] - code[i];
                    prev[i + 1][j] = trellis[i + 1][j].from0;
                } else {
                    if (d[i][trellis[i + 1][j].from1] + code[i] >
                        d[i][trellis[i + 1][j].from0] - code[i]) {
                        d[i + 1][j] = d[i][trellis[i + 1][j].from1] + code[i];
                        prev[i + 1][j] = trellis[i + 1][j].from1;
                    } else {
                        d[i + 1][j] = d[i][trellis[i + 1][j].from0] - code[i];
                        prev[i + 1][j] = trellis[i + 1][j].from0;
                    }
                }
            }
        }
        size_t p = 0;
        uint64_t result = 0;
        for (size_t i = n; i-- > 0;) {
            if (!flip) {
                if (prev[i + 1][p] == trellis[i + 1][p].from1) {
                    result |= (1ULL << i);
                }
            } else {
                if (prev[i + 1][p] == trellis[i + 1][p].from0) {
                    result |= (1ULL << i);
                }
            }
            p = prev[i + 1][p];
        }
        return result;
    }

    void g_msf() {
        size_t starts_fixed = 0;
        for (size_t i = 0; i < k; ++i) {
            replacement[i][i] = 1;
        }

        for (size_t i = 0; i < n && starts_fixed < k; ++i) {
            size_t found_bit = k;
            for (size_t j = starts_fixed; j < k; ++j) {
                if (g[j][i] == 1ULL) {
                    found_bit = j;
                    break;
                }
            }
            if (found_bit == k) {
                continue;
            }
            ++starts_fixed;
            active_start[starts_fixed - 1] = i;
            swap(g[starts_fixed - 1], g[found_bit]);
            swap(replacement[starts_fixed - 1], replacement[found_bit]);
            for (size_t j = starts_fixed; j < k; ++j) {
                if (g[j][i] == 0) {
                    continue;
                }
                for (size_t u = 0; u < n; ++u) {
                    g[j][u] ^= g[found_bit][u];
                }
                for (size_t u = 0; u < k; ++u) {
                    replacement[j][u] ^= replacement[found_bit][u];
                }
            }
        }

        uint64_t ends_fixed = 0;
        for (size_t i = n; i-- > 0 && ends_fixed < k;) {
            size_t found_bit = k;
            for (size_t j = k; j-- > 0;) {
                if (active_end[j] == 0 && g[j][i] == 1ULL) {
                    found_bit = j;
                    break;
                }
            }
            if (found_bit == k) {
                continue;
            }
            ++ends_fixed;
            active_end[found_bit] = i;
            for (size_t j = found_bit; j-- > 0;) {
                if (active_end[j] != 0 || g[j][i] == 0) {
                    continue;
                }
                for (size_t u = 0; u < n; ++u) {
                    g[j][u] ^= g[found_bit][u];
                }
                for (size_t u = 0; u < k; ++u) {
                    replacement[j][u] ^= replacement[found_bit][u];
                }
            }
        }
    }

    uint64_t unfold_by_mask(uint64_t x, uint64_t mask) const {
        uint64_t result = 0;
        size_t p = 0;
        for (size_t i = 0; i < k; ++i) {
            if (((mask >> i) & 1ULL) == 1ULL) {
                result |= (((x >> p) & 1) << i);
                ++p;
            }
        }
        return result;
    }

    uint64_t fold_by_mask(uint64_t x, uint64_t mask) const {
        uint64_t result = 0;
        size_t p = 0;
        for (size_t i = 0; i < k; ++i) {
            if (((mask >> i) & 1ULL) == 1ULL) {
                result |= (((x >> i) & 1) << p);
                ++p;
            }
        }
        return result;
    }

    void build_trellis() {
        uint64_t mask = 0;
        size_t mask_size = 0;
        vector<uint64_t> trellis_mask(n + 1);
        trellis[0].resize(1);
        d[0].resize(1);
        prev[0].resize(1);
        for (size_t i = 1; i < n + 1; ++i) {
            for (size_t j = 0; j < k; ++j) {
                if (active_start[j] == i - 1) {
                    mask |= (1ULL << j);
                    ++mask_size;
                }
                if (active_end[j] == i - 1) {
                    mask &= ~(1ULL << j);
                    --mask_size;
                }
            }
            trellis_mask[i] = mask;
            trellis[i].resize(1ULL << mask_size);
            d[i].resize(1ULL << mask_size);
            prev[i].resize(1ULL << mask_size);
        }

        for (size_t i = 0; i < n; ++i) {
            uint64_t mask_l = trellis_mask[i];
            uint64_t mask_r = trellis_mask[i + 1];
            uint64_t mask_and = mask_l & mask_r;
            for (size_t j1 = 0; j1 < trellis[i].size(); ++j1) {
                for (size_t j2 = 0; j2 < trellis[i + 1].size(); ++j2) {
                    uint64_t values_l = unfold_by_mask(j1, mask_l);
                    uint64_t values_r = unfold_by_mask(j2, mask_r);
                    if (fold_by_mask(values_l, mask_and) == fold_by_mask(values_r, mask_and)) {
                        uint64_t value = 0;
                        uint64_t values_common = values_l | values_r;
                        for (size_t u = 0; u < k; u++) {
                            value ^= ((values_common >> u) & g[u][i]);
                        }
                        if (value == 0) {
                            trellis[i + 1][j2].from0 = j1;
                        } else {
                            trellis[i + 1][j2].from1 = j1;
                        }
                    }
                }
            }
        }
    }
};

uint64_t read_word(size_t k) {
    uint64_t result = 0;
    for (int i = 0; i < k; ++i) {
        uint64_t bit;
        cin >> bit;
        result |= (bit << i);
    }
    return result;
}

void print_word(size_t n, uint64_t word) {
    for (size_t i = 0; i < n; ++i) {
        cout << ((word >> i) & 1) << ' ';
    }
}

int main() {
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    size_t n, k;
    cin >> n >> k;
    vector<vector<uint64_t>> g(k, vector<uint64_t>(n));
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < n; ++j) {
            cin >> g[i][j];
        }
    }
    encoder e(n, k, g);
    coder c(n, k, g);
    for (const auto &x: c.trellis_profile()) {
        cout << x << ' ';
    }
    cout << '\n';

    mt19937 gen;
    uniform_int_distribution<uint64_t> u_dist(0, (1ULL << k) - 1);

    string command;
    while (cin >> command) {
        if (command == "Encode") {
            print_word(n, e.encode(read_word(k)));
        } else if (command == "Decode") {
            vector<double> code(n);
            for (size_t i = 0; i < n; ++i) {
                cin >> code[i];
            }
            print_word(n, c.decode(code, true));
        } else if (command == "Simulate") {
            double noise_level;
            size_t num_of_iterations, max_errors;
            cin >> noise_level >> num_of_iterations >> max_errors;

            normal_distribution<double> n_dist(0, sqrt(0.5 * pow(10, -noise_level / 10) *
                                                       static_cast<double>(n) /
                                                       static_cast<double>(k)));
            size_t total_iterations, total_errors = 0;
            for (total_iterations = 0; total_iterations < num_of_iterations &&
                                       total_errors < max_errors; ++total_iterations) {
                uint64_t word = u_dist(gen);
                uint64_t code = c.encode(word);
                vector<double> noised_code(n);
                for (size_t i = 0; i < n; ++i) {
                    noised_code[i] = static_cast<double>((code >> i) & 1ULL);
                    noised_code[i] = 2 * noised_code[i] - 1;
                    noised_code[i] += n_dist(gen);
                }
                uint64_t received_code = c.decode(noised_code);
                if (code != received_code) {
                    ++total_errors;
                }
            }

            cout << static_cast<double>(total_errors) / static_cast<double>(total_iterations);
        }
        cout << '\n';
    }
}
