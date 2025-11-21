#ifndef SURFACE_TENSION
#define SURFACE_TENSION

template<typename CellList> 
inline bool is_boundary_particle(particles & vd, vect_dist_key_dx a, CellList & NN, double kh)
{
    if (vd.getProp<type>(a) != FLUID) return false;

    Point<2,double> xa = vd.getPos(a);

    std::vector<std::pair<double, double>> covered_arcs;
    
    auto Np = NN.getNNIteratorBox(NN.getCell(xa));
    
    int neighbor_count = 0;
    while (Np.isNext())
    {
        auto b = Np.get();
        
        if (a.getKey() == b) {++Np; continue;};
        
        Point<2,double> xb = vd.getPos(b);
        Point<2,double> dr = xb - xa;
        double r2 = norm2(dr);
        
        if (r2 < kh*kh)
        {
            neighbor_count++;
            double r = sqrt(r2);
            
            double angle = atan2(dr.get(1), dr.get(0));
            if (angle < 0) angle += 2*M_PI;
            
            double theta = acos(r / (2 * H));

            /*double theta = 0;
            if (r <= H) {
                theta = asin(H / r);
            } */

            double alpha1 = angle - theta;
            double alpha2 = angle + theta;
            
            if (alpha1 < 0) alpha1 += 2*M_PI;
            if (alpha2 >= 2*M_PI) alpha2 -= 2*M_PI;
            
            if (alpha1 <= alpha2) {
                covered_arcs.push_back({alpha1, alpha2});
            } else {
                covered_arcs.push_back({alpha1, 2*M_PI});
                covered_arcs.push_back({0, alpha2});
            }
        }
        ++Np;
    }
    
    std::vector<std::pair<double, double>> merged_arcs;
    if (!covered_arcs.empty()) {
        std::sort(covered_arcs.begin(), covered_arcs.end());
        merged_arcs.push_back(covered_arcs[0]);
        
        for (size_t i = 1; i < covered_arcs.size(); i++) {
            auto& last = merged_arcs.back();
            auto& current = covered_arcs[i];
            
            if (current.first <= last.second) {
                last.second = std::max(last.second, current.second);
            } else {
                merged_arcs.push_back(current);
            }
        }
    }
    
    if (merged_arcs.empty()) {
        return true; 
    }
    
    double coverage = 0.0;
    for (const auto& arc : merged_arcs) {
        coverage += (arc.second - arc.first);
    }
    
    return coverage < 2*M_PI - 1e-10;
}

template<typename CellList>
inline Point<2,double> get_new_origin(particles & vd, vect_dist_key_dx a, CellList & NN) {
    double x_sum = 0;
    double y_sum = 0;
    int counter = 0;

    Point<2,double> xa = vd.getPos(a);

    auto Np = NN.getNNIteratorBox(NN.getCell(xa));
    
    while (Np.isNext()) {
        auto b = Np.get();

        if (a.getKey() == b) {
            ++Np;
            continue;
        }

        Point<2,double> xb = vd.getPos(b);

        x_sum += xb.get(0);
        y_sum += xb.get(1);

        counter++;
        ++Np;
    }

    if(counter == 0) {
        return xa;
    }

    double x_new = x_sum / counter;
    double y_new = y_sum / counter;

    return Point<2,double>({x_new, y_new});
}

inline double get_alpha(auto a, const Point<2,double> &new_origin, particles & vd) {
    Point<2,double> xa = vd.getPos(a);

    double rx = xa.get(0) - new_origin.get(0);
    double ry = xa.get(1) - new_origin.get(1);

    double r = sqrt(rx*rx + ry*ry);

    if(rx > 0) {
        return 2 * M_PI - acos(ry / r);
    } else {
        return acos(ry / r);
    }
}

inline Point<2,double> get_local_coord(Point<2,double> xb, const Point<2,double> &new_origin, double alpha) {
    double x_local = (xb.get(0) - new_origin.get(0)) * cos(alpha) + (xb.get(1) - new_origin.get(1)) * sin(alpha);
    double y_local = -(xb.get(0) - new_origin.get(0)) * sin(alpha) + (xb.get(1) - new_origin.get(1)) * cos(alpha);

    return Point<2,double>({x_local, y_local});
}

inline double compute_lagrange(const std::vector<Point<2,double>>& points, double x_eval) {
    double result = 0;
    
    for (size_t j = 0; j < points.size(); j++) {
        double term = points[j].get(1); // y_j
        double xj = points[j].get(0);

        
        for (size_t k = 0; k < points.size(); k++) {
            if (j == k) continue;
            
            double xk = points[k].get(0);
            double denominator = xj - xk;
            
            if (fabs(denominator) < 1e-10) {
                continue;
            }

            term *= (x_eval - xk) / denominator;
        }
        result += term;
    }
    return result;
}

inline double compute_first_deriv_of_lagrange(const std::vector<Point<2,double>>& points, double x_eval) {
    if(points.size() <= 1) {
        return 0;
    }

    double result = 0;
    for (size_t j = 0; j < points.size(); j++) {
        double yj = points[j].get(1);
        double xj = points[j].get(0);
        
        double Lj = 1.0;
        for (size_t k = 0; k < points.size(); k++) {
            if (j == k) continue;
            double xk = points[k].get(0);
            double denominator = xj - xk;
            if (fabs(denominator) < 1e-10) continue;
            Lj *= (x_eval - xk) / denominator;
        }
        
        double sum_reciprocal = 0.0;
        for (size_t m = 0; m < points.size(); m++) {
            if (m == j) continue;
            double xm = points[m].get(0);
            if (fabs(x_eval - xm) < 1e-10) continue; 
            sum_reciprocal += 1.0 / (x_eval - xm);
        }
        
        result += yj * Lj * sum_reciprocal;
    }
    return result;
}

inline double compute_second_deriv_of_lagrange(const std::vector<Point<2,double>>& points, double x_eval) {
    /*if(points.size() <= 2) {
        return 0;
    }

    double result = 0;
    for (size_t j = 0; j < points.size(); j++) {
        double yj = points[j].get(1);
        double xj = points[j].get(0);
        
        double Lj = 1.0;
        for (size_t k = 0; k < points.size(); k++) {
            if (j == k) continue;
            double xk = points[k].get(0);
            double denominator = xj - xk;
            if (fabs(denominator) < 1e-4) continue;
            Lj *= (x_eval - xk) / denominator;
        }
        
        double outer_sum = 0.0;
        
        for (size_t k = 0; k < points.size(); k++) {
            if (k == j) continue;
            
            double inner_product = 1.0;
            
            for (size_t m = 0; m < points.size(); m++) {
                if (m == j || m == k) continue;
                
                double xm = points[m].get(0);
                double diff = x_eval - xm;
                if (fabs(diff) < 1e-4) continue;
                
                inner_product *= 1.0 / diff;
            }
            
            outer_sum += inner_product;
        }
        
        result += yj * Lj * 2.0 * outer_sum;
    }
    return result;*/

    double h = 1e-8; 
    
    double f_plus = compute_lagrange(points, x_eval + h);
    double f_center = compute_lagrange(points, x_eval);
    double f_minus = compute_lagrange(points, x_eval - h);
    
    double second_deriv = (f_plus - 2*f_center + f_minus) / (h * h);
    
    // Сильная регуляризация
    return second_deriv;
}

//refactor me
template<typename CellList>
inline std::vector<Point<2,double>> get_neighbors_on_local_coord(Point<2,double> xa, particles & vd, CellList & NN, Point<2,double> new_origin, double alpha) {
    std::vector<Point<2,double>> local_points;
    
    auto Np_all = NN.getNNIteratorBox(NN.getCell(xa));
    while (Np_all.isNext()) {
        auto neighbor = Np_all.get();

        if (!is_boundary_particle(vd, neighbor, NN, 2 * H)) {
            ++Np_all;
            continue;
        }

        Point<2,double> x_neighbor = vd.getPos(neighbor);
        local_points.push_back(get_local_coord(x_neighbor, new_origin, alpha));
        ++Np_all;
    }

    return local_points;
}

//refactor me
inline std::tuple<double, Point<2,double>> get_curvature_and_normal(double P_deriv, double P_deriv2, double alpha) {
    double kappa = fabs(P_deriv2) / pow(1.0 + P_deriv*P_deriv, 1.5);
    /*if (kappa > 2 / dp) {
        kappa = 2 / dp;
    }*/

    if (kappa > 1000) {
        kappa = 1000;
    }
    
    Point<2,double> n_global;
    if (P_deriv2 < 0) {
        n_global = {
            P_deriv * cos(alpha) + sin(alpha),
            P_deriv * sin(alpha) - cos(alpha)
        };
    } else {
        n_global = {
            -P_deriv * cos(alpha) - sin(alpha),
            -P_deriv * sin(alpha) + cos(alpha)
        };
    }
    
    double norm_global = sqrt(n_global.get(0)*n_global.get(0) + n_global.get(1)*n_global.get(1));
    n_global = {n_global.get(0)/norm_global, n_global.get(1)/norm_global};

    return {kappa, n_global};
}

#endif //SURFACE_TENSION