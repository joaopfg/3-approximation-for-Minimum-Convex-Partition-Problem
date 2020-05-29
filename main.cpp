/*****************************************************************************************
Author: João Fontes Gonçalves

Implementation of the 3-approximation for Minimum Convex Partition Problem

Described in the article:
Christian Knauer and Andreas Spillner: Approximation Algorithms for the Minimum Convex Partition Problem, SWAT 2006, pp. 232-241

It uses this library (https://github.com/nlohmann/json) to parse json files.

It uses this wrapper (https://github.com/lava/matplotlib-cpp) of python's matplotlib library to print the phases of the recursion.

To compile: g++ -o main main.cpp -I/usr/include/python2.7 -lpython2.7
*****************************************************************************************/

#include <bits/stdc++.h>
#include <nlohmann/json.hpp>
#include "matplotlibcpp.h"

using namespace std;
using json = nlohmann::json;

namespace plt = matplotlibcpp;

typedef pair<int,int> ii;

const double eps = 1e-8;

vector<ii> final_sol;
set<ii> seen;

struct pt{
	int ind;
	double x, y;
};

///////////////////////////////////////////////////
struct Triangle{
	pt a, b, c;
};

///////////////////////////////////////////////////

bool Equal(double a, double b){
	return abs(a-b) < eps;
}

///////////////////////////////////////////////////

bool Equal_points(pt a, pt b){
	return Equal(a.x, b.x) && Equal(a.y, b.y);
}

///////////////////////////////////////////////////

bool Equal_edges(pt a, pt b, pt c, pt d){
	return (Equal_points(a,c) && Equal_points(b,d)) || (Equal_points(a,d) && Equal_points(b,c));
}

///////////////////////////////////////////////////

bool cmp(pt a, pt b){
	return a.x < b.x || (Equal(a.x, b.x) && a.y < b.y);
}

///////////////////////////////////////////////////

bool cw(pt a, pt b, pt c){
	return a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y) < 0;
}

///////////////////////////////////////////////////

bool ccw(pt a, pt b, pt c){
	return a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y) > 0;
}

///////////////////////////////////////////////////

double dist(pt a, pt b){
	return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

///////////////////////////////////////////////////

double norm(pt a){
	return sqrt(a.x*a.x + a.y*a.y);
}

///////////////////////////////////////////////////

double dot_product(pt o, pt a, pt b){
	pt oa, ob;
	oa.x = a.x - o.x;
	oa.y = a.y - o.y;
	ob.x = b.x - o.x;
	ob.y = b.y - o.y;
	return oa.x*ob.x + oa.y*ob.y;
}

///////////////////////////////////////////////////

double cos_angle(pt o, pt a, pt b){
	pt oa, ob;
	oa.x = a.x - o.x;
	oa.y = a.y - o.y;
	ob.x = b.x - o.x;
	ob.y = b.y - o.y;
	return dot_product(o, a, b)/(norm(oa)*norm(ob));
}

///////////////////////////////////////////////////

//Test if points p and q are on the same side of line ab
bool same_side(pt a, pt b, pt p, pt q){
	return (cw(a,p,b) && cw(a,q,b)) || (ccw(a,p,b) && ccw(a,q,b));
}

///////////////////////////////////////////////////

//Test if point p is on segment ab
bool point_in_segment(pt p, pt a, pt b){
	return Equal(dist(p,a) + dist(p,b), dist(a,b));
}

///////////////////////////////////////////////////

//Test if point p is on line ab
bool point_in_line(pt p, pt a, pt b){
	return Equal(dist(p,a) + dist(p,b), dist(a,b)) || Equal(abs(dist(p,a) - dist(p,b)), dist(a,b));
}

///////////////////////////////////////////////////

//Test for intersection between two lines ab and cd
//Supposing that the two lines are never coincidents
bool line_intersect(pt a, pt b, pt c, pt d, pt &intersect){
	double m1, h1, m2, h2;

	if(Equal(a.x, b.x) && Equal(c.x, d.x)) return false;
	else if(Equal(a.x, b.x) && !Equal(c.x, d.x)){
		m2 = (d.y - c.y)/(d.x - c.x);
		h2 = d.y - m2*d.x;
		intersect.x = a.x;
		intersect.y = m2*a.x + h2;
		return true;
	}
	else if(!Equal(a.x, b.x) && Equal(c.x, d.x)){
		m1 = (b.y - a.y)/(b.x - a.x);
		h1 = b.y - m1*b.x;
		intersect.x = c.x;
		intersect.y = m1*c.x + h1;
		return true;
	}
	else{
		m1 = (b.y - a.y)/(b.x - a.x);
		h1 = b.y - m1*b.x;
		m2 = (d.y - c.y)/(d.x - c.x);
		h2 = d.y - m2*d.x;

		if(Equal(m1, m2)) return false;
		else{
			intersect.x = (h2 - h1)/(m1 - m2);
			intersect.y = m1*intersect.x + h1;
			return true;
		}
	}
}

///////////////////////////////////////////////////

//Test for intersection between line ab and segment cd
bool line_intersect_segment(pt a, pt b, pt c, pt d, pt &intersect){
	return line_intersect(a, b, c, d, intersect) && point_in_segment(intersect, c, d);
}

///////////////////////////////////////////////////

//Test for intersection between segment ab and segment cd
bool segment_intersect_segment(pt a, pt b, pt c, pt d, pt &intersect){
	return line_intersect(a, b, c, d, intersect) && point_in_segment(intersect, a, b) && point_in_segment(intersect, c, d);
}

///////////////////////////////////////////////////

//Get random integer in the interval [0, up_limit]
int RandUnifInt(int up_limit){
	random_device rd;
    mt19937 gen(rd()); 
    uniform_int_distribution<int> dis(0, up_limit);
    return dis(gen);
}

///////////////////////////////////////////////////

//Search intersection between intern segment ab and segments of the convex hull
//Does not count the possibility of colinear points
void get_intersection(pt a, pt b, pt &r1, pt &q1, pt &r2, pt &q2, vector<pt> &extern_convex_hull, vector<pt> &intern_convex_hull){
	pt intersect;
	pt test = intern_convex_hull[2];
	int cont = 0;
	bool flag1 = false, flag2 = false;

	for(int i=0;i<(int)extern_convex_hull.size() - 1;i++){
		if(line_intersect_segment(a, b, extern_convex_hull[i], extern_convex_hull[i+1], intersect)){
			cont++;

			if(cont == 1){
				if(dist(intersect, a) < dist(intersect, b)){
					if(same_side(a, b, test, extern_convex_hull[i])){
						q1 = extern_convex_hull[i];
						r1 = extern_convex_hull[i+1];
					}
					else{
						r1 = extern_convex_hull[i];
						q1 = extern_convex_hull[i+1];
					}
					flag1 = true;
				}
				else{
					if(same_side(a, b, test, extern_convex_hull[i])){
						q2 = extern_convex_hull[i];
						r2 = extern_convex_hull[i+1];
					}
					else{
						r2 = extern_convex_hull[i];
						q2 = extern_convex_hull[i+1];
					}
					flag2 = true;
				}
			}
			else if(cont == 2){
				if(flag1){
					if(same_side(a, b, test, extern_convex_hull[i])){
						q2 = extern_convex_hull[i];
						r2 = extern_convex_hull[i+1];
					}
					else{
						r2 = extern_convex_hull[i];
						q2 = extern_convex_hull[i+1];
					}
				}
				else{
					if(same_side(a, b, test, extern_convex_hull[i])){
						q1 = extern_convex_hull[i];
						r1 = extern_convex_hull[i+1];
					}
					else{
						r1 = extern_convex_hull[i];
						q1 = extern_convex_hull[i+1];
					}
				}

				break;
			}
		}
	}

	if(cont < 2 && line_intersect_segment(a, b, extern_convex_hull[extern_convex_hull.size() - 1], extern_convex_hull[0], intersect)){
		if(flag1){
			if(same_side(a, b, test, extern_convex_hull[extern_convex_hull.size() - 1])){
				q2 = extern_convex_hull[extern_convex_hull.size() - 1];
				r2 = extern_convex_hull[0];
			}
			else{
				r2 = extern_convex_hull[extern_convex_hull.size() - 1];
				q2 = extern_convex_hull[0];
			}
		}
		else{
			if(same_side(a, b, test, extern_convex_hull[extern_convex_hull.size() - 1])){
				q1 = extern_convex_hull[extern_convex_hull.size() - 1];
				r1 = extern_convex_hull[0];
			}
			else{
				r1 = extern_convex_hull[extern_convex_hull.size() - 1];
				q1 = extern_convex_hull[0];
			}
		}
	}
}

///////////////////////////////////////////////////

//Get convex hull in clock-wise sense
void convex_hull(vector<pt> &a){
	if(a.size() == 0 || a.size() == 1) return;

	sort(a.begin(), a.end(), &cmp);
	pt p1 = a[0], p2 = a.back();
	vector<pt> up, down;

	up.push_back(p1);
	down.push_back(p1);

	for(int i=1;i<(int)a.size();i++){
		if(i == a.size() - 1 || cw(p1,a[i],p2)){
			while(up.size() >= 2 && !cw(up[up.size() - 2], up[up.size() - 1], a[i]))
				up.pop_back();
			up.push_back(a[i]);
		}
		if(i == a.size() - 1 || ccw(p1,a[i],p2)){
			while(down.size() >= 2 && !ccw(down[down.size() - 2], down[down.size() - 1], a[i]))
				down.pop_back();
			down.push_back(a[i]);
		}
	}

	a.clear();
	for(int i=0;i<(int)up.size();i++) a.push_back(up[i]);
	for(int i=down.size() - 2;i>0;i--) a.push_back(down[i]);
}

///////////////////////////////////////////////////

//Plot points using matplotlib wrapper
void plot_points(vector<pt> &points, string &title_plot){
	vector<double> x, y;

	for(int i=0;i<(int)points.size();i++){
		x.push_back(points[i].x);
		y.push_back(points[i].y);
	}

	plt::title(title_plot);
    plt::scatter(x, y, 10.0);
    plt::show();
}

///////////////////////////////////////////////////

//Get triangle area
double triangle_area(Triangle &t){
	return abs((t.b.x - t.a.x)*(t.c.y - t.a.y) - (t.b.y - t.a.y)*(t.c.x - t.a.x))/2.0;
}

///////////////////////////////////////////////////

//Test if a point is inside a triangle
bool inside_triangle(pt p, Triangle &t){
	Triangle t1, t2, t3;
	t1.a = t.a;
	t1.b = t.b;
	t1.c = p;
	t2.a = t.a;
	t2.b = t.c;
	t2.c = p;
	t3.a = t.b;
	t3.b = t.c;
	t3.c = p;

	return Equal(triangle_area(t1) + triangle_area(t2) + triangle_area(t3), triangle_area(t));
}

///////////////////////////////////////////////////

//Get the convex hull area (useful to test if a point is inside the convex hull)
double convex_hull_area(vector<pt> convex_hull){
	double area = 0.0;

	for(int i=1;i<(int)convex_hull.size()-1;i++){
		Triangle t;
		t.a = convex_hull[0];
		t.b = convex_hull[i];
		t.c = convex_hull[i+1];
		area += triangle_area(t);
	}

	return area; 
}

///////////////////////////////////////////////////

//Test if a point is inside the convex hull
bool inside_convex_hull(pt p, vector<pt> &convex_hull){
	double area_point = 0.0;
	Triangle t;

	for(int i=0;i<(int)convex_hull.size() - 1;i++){	
		t.a = p;
		t.b = convex_hull[i];
		t.c = convex_hull[i+1];
		area_point += triangle_area(t);
	}

	t.a = p;
	t.b = convex_hull[convex_hull.size() - 1];
	t.c = convex_hull[0];
	area_point += triangle_area(t);

	return Equal(area_point, convex_hull_area(convex_hull));
}

///////////////////////////////////////////////////

//Get a random triangle with the vertices of the convex hull
Triangle random_triangle(vector<pt> &convex_hull){
	int sz = convex_hull.size() - 1;

	int inda = RandUnifInt(sz), indb, indc;

	while(true){
		indb = RandUnifInt(sz);

		if(indb != inda) break;
	}

	while(true){
		indc = RandUnifInt(sz);

		if(indc != inda && indc != indb) break;
	}

	Triangle t;
	t.a = convex_hull[inda];
	t.b = convex_hull[indb];
	t.c = convex_hull[indc];
	return t;
}

///////////////////////////////////////////////////

//Solve case in which there is zero points in the interior
//Suppose convex hull points given in clockwise order
vector<ii> solve_case0(vector<pt> &points){
	vector<ii> inds;

	for(int i=0;i<(int)points.size()-1;i++)
		inds.push_back({points[i].ind, points[i+1].ind});

	inds.push_back({points[points.size()-1].ind, points[0].ind});

	return inds;
}

///////////////////////////////////////////////////

//Solve case in which there is only one point in the interior
vector<ii> solve_case1(vector<pt> &points){
	vector<pt> points_copy;
	pt p;

	for(int i=0;i<(int)points.size();i++) points_copy.push_back(points[i]);

	convex_hull(points);

	set<int> convex_hull_inds;

	for(int i=0;i<(int)points.size();i++) convex_hull_inds.insert(points[i].ind);

	for(int i=0;i<(int)points_copy.size();i++){
		if(convex_hull_inds.count(points_copy[i].ind) == 0){
			p = points_copy[i];
			break;
		} 
	}

	Triangle t;
	vector<ii> edges_hull = solve_case0(points);

	while(true){
		t = random_triangle(points);
		if(inside_triangle(p, t)) break; 
	}

	edges_hull.push_back({p.ind, t.a.ind});
	edges_hull.push_back({p.ind, t.b.ind});
	edges_hull.push_back({p.ind, t.c.ind});

	return edges_hull;
}

///////////////////////////////////////////////////

//Given segment qr exterior to the convex hull
//Get the point p in the convex hull such that the angle between qp and qr is the minimal one
pt get_tangent_to_convex_hull(pt q, pt r, vector<pt> &convex_hull){
	double cos_max = -2.0;
	int ind_max;

	for(int i=0;i<(int)convex_hull.size();i++){
		if(cos_angle(q, r, convex_hull[i]) > cos_max){
			cos_max = cos_angle(q, r, convex_hull[i]);
			ind_max =  i;
		}
	}

	return convex_hull[ind_max];
}

///////////////////////////////////////////////////

//Recursion to solve the general case
//Suppose that there is no three points in a same line
void solve(vector<pt> &points){
	vector<pt> points_copy, points_inside, points_inside_copy;

	for(int i=0;i<(int)points.size();i++) points_copy.push_back(points[i]);

	convex_hull(points);

	cout << "Convex hull: ";
	for(int i=0;i<(int)points.size();i++){
		cout << points[i].ind << " ";
	}
	cout << endl;

	string title_plot = "Next phase convex hull";
	plot_points(points, title_plot);

	set<int> convex_hull_inds;

	for(int i=0;i<(int)points.size();i++) convex_hull_inds.insert(points[i].ind);

	for(int i=0;i<(int)points_copy.size();i++){
		if(convex_hull_inds.count(points_copy[i].ind) == 0) points_inside.push_back(points_copy[i]); 
	}

	for(int i=0;i<(int)points_inside.size();i++) points_inside_copy.push_back(points_inside[i]);

	convex_hull(points_inside);

	cout << "Convex hull inside: ";
	for(int i=0;i<(int)points_inside.size();i++)
		cout << points_inside[i].ind << " ";
	cout << endl;

	if(points_inside_copy.size() == 0){
		cout << "entry0" << endl;
		vector<ii> edges0 = solve_case0(points);

		for(int i=0;i<(int)edges0.size();i++){
			if(seen.count(edges0[i]) == 0){
				seen.insert(edges0[i]);
				final_sol.push_back(edges0[i]);
			}
		} 

		return;
	}
	else if(points_inside_copy.size() == 1){
		cout << "entry1" << endl;
		vector<ii> edges1 = solve_case1(points_copy);

		for(int i=0;i<(int)edges1.size();i++){
			if(seen.count(edges1[i]) == 0){
				seen.insert(edges1[i]);
				final_sol.push_back(edges1[i]);
			}
		} 

		return;
	}
	else{
		int ind = 0, broken_pos = -1;
		pt r1, q1, r2, q2;
		vector<int> necklace;
		set<int> necklace_set;
		bool zero = false;

		get_intersection(points_inside[ind], points_inside[ind+1], r1, q1, r2, q2, points, points_inside);

		pt t1 = get_tangent_to_convex_hull(q1, r1, points_inside);
		pt t2 = get_tangent_to_convex_hull(q2, r2, points_inside);

		if(seen.count({q1.ind, r1.ind}) == 0){
			seen.insert({q1.ind, r1.ind});
			final_sol.push_back({q1.ind, r1.ind});
		}

		if(seen.count({q2.ind, r2.ind}) == 0){
			seen.insert({q2.ind, r2.ind});
			final_sol.push_back({q2.ind, r2.ind});
		}

		if(seen.count({r1.ind, t1.ind}) == 0){
			seen.insert({r1.ind, t1.ind});
			final_sol.push_back({r1.ind, t1.ind});
		}
		
		if(seen.count({r2.ind, t2.ind}) == 0){
			seen.insert({r2.ind, t2.ind});
			final_sol.push_back({r2.ind, t2.ind});
		}
		
		if(seen.count({r1.ind, points_inside[ind].ind}) == 0){
			seen.insert({r1.ind, points_inside[ind].ind});
			final_sol.push_back({r1.ind, points_inside[ind].ind});
		}
		
		if(seen.count({r2.ind, points_inside[ind+1].ind}) == 0){
			seen.insert({r2.ind, points_inside[ind+1].ind});
			final_sol.push_back({r2.ind, points_inside[ind+1].ind});
		}

		for(int i=0;i<(int)points.size();i++){
			//If point is between r1 and r2
			if((!Equal_points(points[i], q1) && !Equal_points(points[i], q2) &&
			 cos_angle(points[i], r1, r2) < cos_angle(points[i], q1, q2)) ||
			 Equal_points(points[i], r1) || Equal_points(points[i], r2)){
				necklace.push_back(i);
				necklace_set.insert(i);

				if(i == 0) zero = true;
			}
		}

		cout << "Necklace: ";
		for(int i=0;i<(int)necklace.size();i++){
			cout << points[necklace[i]].ind << " ";
		}
		cout << endl;

		if(zero){
			for(int i=1;i<(int)necklace.size();i++){
				if(necklace[i] != necklace[i-1] + 1){
					broken_pos = i;
					break;
				}
			}

			if(broken_pos != -1){
				for(int i=0;i<broken_pos-1;i++){
					if(seen.count({points[necklace[i]].ind, points[necklace[i+1]].ind}) == 0){
						seen.insert({points[necklace[i]].ind, points[necklace[i+1]].ind});
						final_sol.push_back({points[necklace[i]].ind, points[necklace[i+1]].ind});
					}
				}

				for(int i=broken_pos;i<(int)necklace.size() - 1;i++){
					if(seen.count({points[necklace[i]].ind, points[necklace[i+1]].ind}) == 0){
						seen.insert({points[necklace[i]].ind, points[necklace[i+1]].ind});
						final_sol.push_back({points[necklace[i]].ind, points[necklace[i+1]].ind});
					}
				}

				if(seen.count({points[necklace[broken_pos-1]].ind, points[necklace[broken_pos]].ind}) == 0){
					seen.insert({points[necklace[broken_pos-1]].ind, points[necklace[broken_pos]].ind});
					final_sol.push_back({points[necklace[broken_pos-1]].ind, points[necklace[broken_pos]].ind});
				}
			}
		}
		
		if(broken_pos == -1){
			for(int i=0;i<(int)necklace.size() - 1;i++){
				if(seen.count({points[necklace[i]].ind, points[necklace[i+1]].ind}) == 0){
					seen.insert({points[necklace[i]].ind, points[necklace[i+1]].ind});
					final_sol.push_back({points[necklace[i]].ind, points[necklace[i+1]].ind});
				}	
			}
		}

		vector<pt> new_convex_hull;

		for(int i=0;i<(int)points.size();i++){
			if(!Equal_points(points[i], r1) && !Equal_points(points[i], r2) && necklace_set.count(i) == 0) 
				new_convex_hull.push_back(points[i]);
		}

		new_convex_hull.push_back(t1);
		new_convex_hull.push_back(t2);
		new_convex_hull.push_back(points_inside[ind]);
		new_convex_hull.push_back(points_inside[ind+1]);

		convex_hull(new_convex_hull);

		vector<pt> new_points;

		for(int i=0;i<(int)points_copy.size();i++){
			if(inside_convex_hull(points_copy[i], new_convex_hull))
				new_points.push_back(points_copy[i]);
		}

		solve(new_points);

		return;
	}
}

///////////////////////////////////////////////////

//Calculate score of the solution
double get_score(vector<pt> &points){
	vector<pt> convex_hull_points;

	for(int i=0;i<(int)points.size();i++) 
		convex_hull_points.push_back(points[i]);

	convex_hull(convex_hull_points);

	int n = points.size();
	int c = convex_hull_points.size();
	int m = final_sol.size();
	int s = 3*(n-1) - c - m;
	double score = (1.0*s)/(1.0*(3*(n-1) - c));

	return score;
}

///////////////////////////////////////////////////

//Test if point set has three collinear points
bool three_collinear(vector<pt> &points){
	for(int i=0;i<(int)points.size();i++){
		for(int j=i+1;j<(int)points.size();j++){
			for(int k=j+1;k<(int)points.size();k++){
				if(point_in_line(points[k], points[i], points[j])) return true;
			}
		}
	}

	return false;
}

///////////////////////////////////////////////////

//Verify if the solution is ok
bool valid_solution(vector<pt> &points){
	int n = points.size();
	pt intersect;
	//No edge is loop that connects a vertex to itself
	for(int i=0;i<(int)final_sol.size();i++){
		if(final_sol[i].first == final_sol[i].second) return false;
	}

	cout << "Passed test 1" << endl;

	//Each vertex is incident to at least two edges
	vector<int> cont(n);
	for(int i=0;i<(int)final_sol.size();i++){
		cont[final_sol[i].first]++;
		cont[final_sol[i].second]++;
	}

	for(int i=0;i<n;i++){
		if(cont[i] < 2){
			cout << i << endl;
			return false;
		} 
	}

	cout << "Passed test 2" << endl;

	//No edge intersects another edge in the interior
	for(int i=0;i<(int)final_sol.size();i++){
		for(int j=i+1;j<(int)final_sol.size();j++){
			if(segment_intersect_segment(points[final_sol[i].first], points[final_sol[i].second], points[final_sol[j].first], points[final_sol[j].second], intersect) &&
				!Equal_points(points[final_sol[i].first], points[final_sol[j].first]) &&
				!Equal_points(points[final_sol[i].first], points[final_sol[j].second]) &&
				!Equal_points(points[final_sol[i].second], points[final_sol[j].first]) &&
				!Equal_points(points[final_sol[i].second], points[final_sol[j].second])) return false;
		}
	}

	cout << "Passed test 3" << endl;

	//No point lies in the interior of an edge
	for(int i=0;i<(int)points.size();i++){
		for(int j=0;j<(int)final_sol.size();j++){
			if(point_in_segment(points[i], points[final_sol[j].first], points[final_sol[j].second]) &&
				!Equal_points(points[i], points[final_sol[j].first]) &&
				!Equal_points(points[i], points[final_sol[j].second])) return false;
		}
	}

	cout << "Passed test 4" << endl;
	//The boundary of the convex hull of the points is covered (exactly) by a subset of output edges
	vector<pt> convex_hull_points;

	for(int i=0;i<(int)points.size();i++) 
		convex_hull_points.push_back(points[i]);

	convex_hull(convex_hull_points);

	bool good;

	for(int i=0;i<(int)convex_hull_points.size() - 1;i++){
		good = false;

		for(int j=0;j<(int)final_sol.size();j++){
			if(Equal_edges(convex_hull_points[i], convex_hull_points[i+1], points[final_sol[j].first], points[final_sol[j].second])){
				good = true;
				break;
			}
		}

		if(!good) return false;
	}

	good = false;

	for(int j=0;j<(int)final_sol.size();j++){
		if(Equal_edges(convex_hull_points[convex_hull_points.size() - 1], convex_hull_points[0], points[final_sol[j].first], points[final_sol[j].second])){
			good = true;
			break;
		}
	}

	if(!good) return false;

	cout << "Passed test 5" << endl;

	return true;
}

int main()
{	
	string s;
    ifstream file_in("images/euro-night-0000025.instance.json");
    while(file_in.is_open()){
    	getline(file_in, s);
    	file_in.close();
    }

    //Parse and serialize JSON
    json j_complete = json::parse(s);
    vector<pt> points, points_copy;

    for(int i=0;i<(int)j_complete["points"].size();i++){
    	pt point;
    	point.ind = i;
    	point.x = j_complete["points"][i]["x"];
    	point.y = j_complete["points"][i]["y"];
    	points.push_back(point);
    }

    for(int i=0;i<(int)points.size();i++)
    	points_copy.push_back(points[i]);

    if(three_collinear(points_copy)){
    	cout << "Three points are collinear" << endl;
    	cout << "Not possible to compute!" << endl;
    	return 0;
    } 
    else cout << "No three points are collinear" << endl;

    string text = "Points";
    plot_points(points_copy, text);

    solve(points);

    cout << "The final solution is: " << endl;
    for(int i=0;i<(int)final_sol.size();i++) 
    	cout << "(" << final_sol[i].first << ", " << final_sol[i].second << ")" << endl;

    if(valid_solution(points_copy)) cout << "The solution is valid!" << endl;

    cout << "Score of the solution is: " << get_score(points_copy) << endl;

    return 0;
}