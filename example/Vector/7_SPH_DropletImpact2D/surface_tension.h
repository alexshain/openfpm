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
    

    bool is_boundary = coverage < 2*M_PI - 1e-10;
    vd.template getProp<isBoundary>(a) = is_boundary ? 1 : 0;
    return is_boundary;
}

inline bool weighted_quadratic_fit(const std::vector<Point<2,double>>& pts,
                                  double x_eval,
                                  double &P, double &P1, double &P2) {
    size_t n = pts.size();
    if (n < 3) return false;

    double A00 = 0, A01 = 0, A02 = 0, A11 = 0, A12 = 0, A22 = 0;
    double b0 = 0, b1 = 0, b2 = 0;

    for (size_t i = 0; i < n; ++i) {
        double x = pts[i].get(0);
        double y = pts[i].get(1);
        double dx = x - x_eval;
        double r = fabs(dx);
        
        double w = Wab(r);
        if (w < 1e-12) continue;

        double X0 = 1.0, X1 = x, X2 = x*x;
        
        A00 += w * X0 * X0;
        A01 += w * X0 * X1;
        A02 += w * X0 * X2;
        A11 += w * X1 * X1;
        A12 += w * X1 * X2;
        A22 += w * X2 * X2;

        b0 += w * X0 * y;
        b1 += w * X1 * y;
        b2 += w * X2 * y;
    }

    double det = A00*(A11*A22 - A12*A12) - A01*(A01*A22 - A12*A02) + A02*(A01*A12 - A11*A02);
    if (fabs(det) < 1e-14) return false;

    double det0 = b0*(A11*A22 - A12*A12) - A01*(b1*A22 - A12*b2) + A02*(b1*A12 - A11*b2);
    double det1 = A00*(b1*A22 - A12*b2) - b0*(A01*A22 - A12*A02) + A02*(A01*b2 - b1*A02);
    double det2 = A00*(A11*b2 - A12*b1) - A01*(A01*b2 - A12*b0) + A02*(A01*b1 - A11*b0);

    double a0 = det0 / det;
    double a1 = det1 / det;
    double a2 = det2 / det;

    P = a0 + a1 * x_eval + a2 * x_eval * x_eval;
    P1 = a1 + 2.0 * a2 * x_eval;
    P2 = 2.0 * a2;
    
    return true;
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

inline double get_alpha(const Point<2,double>& xa, const Point<2,double> &new_origin, particles & vd) {
    double rx = xa.get(0) - new_origin.get(0);
    double ry = xa.get(1) - new_origin.get(1);
    
    double theta = atan2(ry, rx);
    
    double alpha = theta - M_PI_2;
    
    if (alpha < 0) alpha += 2*M_PI;
    return alpha;
}

inline Point<2,double> get_local_coord(Point<2,double> xb, const Point<2,double> &new_origin, double alpha) {
    double x_local = (xb.get(0) - new_origin.get(0)) * cos(alpha) + (xb.get(1) - new_origin.get(1)) * sin(alpha);
    double y_local = -(xb.get(0) - new_origin.get(0)) * sin(alpha) + (xb.get(1) - new_origin.get(1)) * cos(alpha);

    return Point<2,double>({x_local, y_local});
}

//refactor me
template<typename CellList>
inline std::vector<Point<2,double>> get_neighbors_on_local_coord(Point<2,double> xa, particles & vd, CellList & NN, Point<2,double> new_origin, double alpha) {
    std::vector<Point<2,double>> local_points;
    
    auto Np_all = NN.getNNIteratorBox(NN.getCell(xa));

    while (Np_all.isNext()) {
        auto neighbor = Np_all.get();

        if (vd.template getProp<isBoundary>(neighbor) == 0) {
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
    if (kappa > 2 / dp) {
        kappa = 2 / dp;
    }

    /*if (kappa > 1000) {
        kappa = 1000;
    }*/
    
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