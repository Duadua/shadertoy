vec4 color_arr[100];
int color_num = 0;
void add_color(vec4 c) { if(color_num < 100) color_arr[color_num++] = c; }
void cac_color(out vec4 c) {
	c = vec4(0.0);
	for(int i = 0; i < color_num; ++i) {
		c = mix(c, color_arr[i], color_arr[i].a);
	}
}

vec2 norm_uv(vec2 uv) { return (2.0 * uv - iResolution.xy) / iResolution.y; }

float cross2(vec2 a, vec2 b) { return a.x*b.y - a.y*b.x; }	// 二维叉积

// ============================================================
// 2d_sdf -- singed distance field
float sd_circle(vec2 p, float r) { return length(p) - r; }

float sd_line(vec2 p, vec2 a, vec2 b, float w) {
	vec2 pa = p - a, ba = b - a;
	return abs(cross2(pa, ba)) / length(ba) - w / 2.0;
}
float sd_seg_line(vec2 p, vec2 a, vec2 b, float w) {
	vec2 pa = p - a, ba = b - a;
	float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
	return length(pa - ba * h) - w / 2.0;
}
float sd_bezier_line(vec2 pos, vec2 A, vec2 B, vec2 C, float w) {			// bezier 曲线
	// A C 为端点, B 为 bezier 点

	vec2 a = B - A;
    vec2 b = A - 2.0*B + C;
    vec2 c = a * 2.0;
    vec2 d = A - pos;

    float kk = 1.0 / dot(b, b);
    float kx = kk * dot(a, b);
    float ky = kk * (2.0*dot(a, a) + dot(d, b)) / 3.0;
    float kz = kk * dot(d, a);      

    float res = 0.0;

    float p = ky - kx*kx;
    float p3 = p*p*p;
    float q = kx*(2.0*kx*kx - 3.0*ky) + kz;
    float h = q*q + 4.0*p3;

    if(h >= 0.0) { 
        h = sqrt(h);
        vec2 x = (vec2(h, -h) - q) / 2.0;
        vec2 uv = sign(x)*pow(abs(x), vec2(1.0 / 3.0));
        float t = uv.x + uv.y - kx;
        t = clamp(t, 0.0, 1.0);

        vec2 qos = d + (c + b*t)*t;
        res = dot(qos,qos);
    }
    else {
        float z = sqrt(-p);
        float v = acos(q / (p*z*2.0)) / 3.0;
        float m = cos(v);
        float n = sin(v)*1.732050808;
        vec3 t = vec3(m + m, -n - m, n - m) * z - kx;
        t = clamp(t, 0.0, 1.0);

        vec2 qos = d + (c + b*t.x)*t.x;
        res = dot(qos, qos);

        qos = d + (c + b*t.y)*t.y;
        res = min(res, dot(qos, qos));

        qos = d + (c + b*t.z)*t.z;
        res = min(res, dot(qos, qos));
    }
    
    return sqrt(res) - w / 2.0;
}

float sd_triangle(vec2 p, vec2 p0, vec2 p1, vec2 p2) {
	vec2 e0 = p1 - p0, e1 = p2 - p1, e2 = p0 - p2;
	vec2 v0 = p  - p0, v1 = p  - p1, v2 = p  - p2;

	vec2 pq0 = v0 - e0*clamp(dot(v0, e0) / dot(e0, e0), 0.0, 1.0);
	vec2 pq1 = v1 - e1*clamp(dot(v1, e1) / dot(e1, e1), 0.0, 1.0);
	vec2 pq2 = v2 - e2*clamp(dot(v2, e2) / dot(e2, e2), 0.0, 1.0);

	float s = sign(cross2(e0, e2));
	vec2 d = min(min(vec2(dot(pq0, pq0), s*(cross2(v0, e0))),
					 vec2(dot(pq1, pq1), s*(cross2(v1, e1)))),
					 vec2(dot(pq2, pq2), s*(cross2(v2, e2))));
	return -sqrt(d.x)*sign(d.y);
}
float sd_isos_triangle(vec2 p, vec2 q) {						// 等腰
	p.x = abs(p.x);
	q.y = -q.y;
    p.y += q.y;

	vec2 a = p - q*clamp(dot(p, q) / dot(q, q), 0.0, 1.0);
	vec2 b = p - q*vec2(clamp(p.x / q.x, 0.0, 1.0), 1.0);
	float s = -sign(q.y);
	vec2 d = min(vec2(dot(a, a), s*(cross2(p, q))), vec2(dot(b, b), s*(p.y - q.y)));
	return -sqrt(d.x) * sign(d.y);
}
float sd_equi_triangle(vec2 p, float r) {						// 等边
	const float k = sqrt(3.0);

	p.x = abs(p.x) - r;
	if(p.x + k*p.y > 0.0) p = vec2(p.x - k*p.y, -k*p.x - p.y) / 2.0;
	p.x -= clamp(p.x, -2.0, 0.0);
	return -length(p)*sign(p.y);
}

float sd_rect(vec2 p, vec2 size) {
	vec2 d = abs(p) - size;										// abs(vec2) -- 每个分量取绝对值
	return length(max(d, vec2(0.0))) + min(max(d.x, d.y), 0.0);	// max(vec2, vec2) -- 每个分量取最大值 
} // 取差向量的 正分量向量长度 + 最大的负分量

float sd_rhombus(vec2 p, vec2 b) {							// 菱形
	vec2 q = abs(p);
	
	vec2 c = 0.5*b;
	vec2 e = q - c;
	c.x = -c.x;
	float h = clamp(dot(e, c) / dot(c, c), -1.0, 1.0);
	float d = length(e - c*h);
	return d*sign(cross2(e, c));
}

float sd_ellipse(vec2 p, vec2 ab) {					// 椭圆
    p = abs(p); if( p.x > p.y ) { p = p.yx; ab = ab.yx; }
    float l = ab.y*ab.y - ab.x*ab.x;
	
    float m = ab.x*p.x/l;      float m2 = m*m; 
    float n = ab.y*p.y/l;      float n2 = n*n; 
    float c = (m2 + n2 - 1.0) / 3.0; float c3 = c*c*c;
	
    float q = c3 + m2*n2*2.0;
    float d = c3 + m2*n2;
    float g = m + m*n2;

    float co;
    if(d < 0.0) {
        float h = acos(q / c3) / 3.0;
        float s = cos(h);
        float t = sin(h)*sqrt(3.0);
        float rx = sqrt(-c*(s + t + 2.0) + m2);
        float ry = sqrt(-c*(s - t + 2.0) + m2);
        co = (ry + sign(l)*rx + abs(g)/(rx*ry) - m) / 2.0;
    }
    else {
        float h = 2.0*m*n*sqrt(d);
        float s = sign(q + h)*pow(abs(q + h), 1.0 / 3.0);
        float u = sign(q - h)*pow(abs(q - h), 1.0 / 3.0);
        float rx = -s - u - c*4.0 + 2.0*m2;
        float ry = (s - u)*sqrt(3.0);
        float rm = sqrt(rx*rx + ry*ry);
        co = (ry / sqrt(rm - rx) + 2.0*g / rm - m) / 2.0;
    }

    vec2 r = ab*vec2(co, sqrt(1.0 - co*co));
    return length(r - p) * sign(p.y - r.y);
}

float sd_isos_trapezoid(vec2 p, float r1, float r2, float h) {	// 等腰梯形
	vec2 k1 = vec2(r2, h);
    vec2 k2 = vec2(r2-r1, 2.0*h);

    p.x = abs(p.x);
    vec2 ca = vec2(p.x - min(p.x, (p.y < 0.0) ? r1 : r2), abs(p.y) - h);
    vec2 cb = p - k1 + k2*clamp(dot(k1 - p, k2) / dot(k2, k2), 0.0, 1.0 );
    
    float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
    
    return s*sqrt(min(dot(ca, ca), dot(cb, cb)));
}

float sd_uneven_capsule(vec2 p, float r1, float r2, float h) {	// 不均匀胶囊体
    p.x = abs(p.x);
    
    float b = (r1 - r2) / h;
    float a = sqrt(1.0 - b*b);
    float k = dot(p, vec2(-b, a));
    
    if(k < 0.0) return length(p) - r1;
    if(k > a*h) return length(p - vec2(0.0, h)) - r2;
        
    return dot(p, vec2(a,b)) - r1;
}

float sd_pentagon(vec2 p, float r) {						// 正五边形
    const vec3 k = vec3(0.809016994,0.587785252,0.726542528);
    p.x = abs(p.x);
    p -= 2.0*min(dot(vec2(-k.x, k.y), p), 0.0)*vec2(-k.x, k.y);
    p -= 2.0*min(dot(vec2( k.x, k.y), p), 0.0)*vec2( k.x, k.y);
    return length(p - vec2(clamp(p.x, -r*k.z, r*k.z), r))*sign(p.y - r);
}
float sd_octogon(in vec2 p, in float r) {					// 正八边形
  	const vec3 k = vec3(-0.9238795325, 0.3826834323, 0.4142135623 );
  	p = abs(p);
  	p -= 2.0*min(dot(vec2( k.x, k.y), p), 0.0)*vec2( k.x, k.y);
  	p -= 2.0*min(dot(vec2(-k.x, k.y), p), 0.0)*vec2(-k.x, k.y);
  	return length(p - vec2(clamp(p.x, -k.z*r, k.z*r), r))*sign(p.y - r);
}

float sd_vesica(vec2 p, float r, float d) {					// vesica -- r > d
	p = abs(p);

    float b = sqrt(r*r-d*d); 
    return ((p.y - b)*d > p.x*b) ? length(p - vec2(0.0, b)) : length(p - vec2(-d, 0.0)) - r;
}

float sd_cross(vec2 p, vec2 b, float r) {					// cross
	p = abs(p); p = (p.y > p.x) ? p.yx: p.xy;

	vec2 q = p - b;
	float k = max(q.y, q.x);
	vec2 w = (k > 0.0) ? q : vec2(b.y - p.x, -k);
	return sign(k)*length(max(w, 0.0)) + r;
}

// ============================================================
float sd_round(float sd, float r) { return sd - r; }			// 圆角
float sd_annular(float sd, float r) {return abs(sd) - r; }		// 圆角 + 空心

vec4 draw_2d(float sd, vec3 color, float att) {
	float t = smoothstep(0.0, att, sd);
	return vec4(color, 1.0 - t);
}

// ============================================================

vec4 draw_triangle(vec2 p, vec2 o, vec2 a, vec2 b, vec2 c, vec4 color) {
	vec2 ab = b - a, ac = c - a, bc = c - b;
	
	p = p - o;
	vec2 pa = a - p, pb = b - p, pc = c - p;
	float r1 = pa.x * pb.y - pb.x * pa.y;
	float r2 = pb.x * pc.y - pc.x * pb.y;
	float r3 = pc.x * pa.y - pa.x * pc.y;
	if(r1 > 0.0 && r2 > 0.0 && r3 > 0.0) return color;
	if(r1 < 0.0 && r2 < 0.0 && r3 < 0.0) return color;
	return vec4(0.0);
}

vec4 draw_rect(vec2 p, vec2 o, vec2 size, vec4 color) {
	vec4 r = vec4(0.0);
	vec2 a = vec2(0.0), b = vec2(size.x, 0.0), c = vec2(0.0, size.y);
	r =	draw_triangle(p, o, a, b, c, color);
	if(length(r) > 0.0) return r;
	r = draw_triangle(p, o, b, c, size, color);
	return r;

}

void mainImage(out vec4 frag_color, in vec2 frag_coord) {
	
	vec2 uv = norm_uv(frag_coord);

	vec2 ro = norm_uv(vec2(iMouse.xy));

	// draw point 
	vec2 o = vec2(0.0, 0.0);
	vec2 o1 = vec2(0.0, 1.0);
	vec2 o2 = vec2(1.0, 0.0);
	vec4 c_circle = draw_2d(sd_circle(uv - (ro + o), 0.02), vec3(0.4,0.2,0.6), 0.01);
	vec4 c_circle1 = draw_2d(sd_circle(uv - (ro + o1), 0.02), vec3(0.4,0.2,0.6), 0.01);
	vec4 c_circle2 = draw_2d(sd_circle(uv - (ro + o2), 0.02), vec3(0.4,0.2,0.6), 0.01);

	// draw circle
	vec2 o3 = vec2(1.0, 0.5);
	vec4 c_circle3 = draw_2d(sd_circle(uv - (ro + o3), 0.2), vec3(0.4,0.0,0.2), 0.01);
	// draw ellipse
	vec2 o_ellipse = vec2(1.5, 0.5);
	vec4 c_ellipse = draw_2d(sd_ellipse(uv - (ro + o_ellipse), vec2(0.2, 0.1)), vec3(0.4,0.0,0.2), 0.01);

	// draw line
	vec4 c_line1 = draw_2d(sd_line(uv, ro + o, ro + o1, 0.01), vec3(0.5), 0.01);
	vec4 c_line2 = draw_2d(sd_line(uv, ro + o, ro + o2, 0.01), vec3(0.5), 0.01);
	// draw segment line
	vec4 c_seg_line1 = draw_2d(sd_seg_line(uv, ro + o, ro + o1, 0.01), vec3(0.5, 1.0, 0.2), 0.01);
	vec4 c_seg_line2 = draw_2d(sd_seg_line(uv, ro + o, ro + o2, 0.01), vec3(0.5, 1.0, 0.2), 0.01);
	vec4 c_seg_line3 = draw_2d(sd_seg_line(uv, ro + o1, ro + o2, 0.01), vec3(0.5, 1.0, 0.2), 0.01);
	// draw bezier line
	vec4 c_bezier_line = draw_2d(sd_bezier_line(uv, ro + o1, ro + o, ro + o2, 0.01), vec3(0.4,0.0,0.2), 0.01);
	
	// draw_triangle
	vec2 o_equi_triangle = vec2(-1.5, -0.9);
	vec2 o_isos_triangle = vec2(-1.3, -0.9);
	vec2 o_triangle = vec2(-1.0, -0.9);
	vec2 A = vec2(-0.3, 0.5);
	vec2 B = vec2(0.5, 0.3);
	vec2 C = vec2(-0.0, -0.0);
	vec4 c_triangle1 = draw_2d(sd_equi_triangle(uv - (ro + o_equi_triangle), 0.1), vec3(0.0, 0.5, 0.5), 0.01);
	vec4 c_triangle2 = draw_2d(sd_isos_triangle(uv - (ro + o_isos_triangle), vec2(0.1, 0.4)), vec3(0.0, 0.5, 0.5), 0.01);
	vec4 c_triangle3 = draw_2d(sd_triangle(uv - (ro + o_triangle), A, B, C), vec3(0.0, 0.5, 0.5), 0.01);

	// draw rect 
	vec2 o_rect = vec2(-1.0, 0.6);
	vec4 c_rect = draw_2d(sd_rect(uv - (ro + o_rect), vec2(0.3, 0.2)), vec3(1.0, 0.0, 0.0), 0.01);
	// draw rhombus
	vec2 o_rhombus = vec2(-1.5, 0.6);
	vec4 c_rhombus = draw_2d(sd_rhombus(uv - (ro + o_rhombus), vec2(0.1, 0.2)), vec3(0.0, 0.0, 1.0), 0.01);

	// draw trapezoid
	vec2 o_trapezoid = vec2(0.2, -0.2);
	vec4 c_trapezoid = draw_2d(sd_isos_trapezoid(uv - (ro + o_trapezoid), 0.1, 0.05, 0.1), vec3(0.5, 0.5, 0.0), 0.01);
	// draw uneven capsule
	vec2 o_uneven_capsule = vec2(0.5, -0.2);
	vec4 c_uneven_capsule = draw_2d(sd_uneven_capsule(uv - (ro + o_uneven_capsule), 0.1, 0.05, 0.1), vec3(0.5, 0.5, 0.0), 0.01);

	// draw pentagon
	vec2 o_pentagon = vec2(0.8, -0.2);
	vec4 c_pentagon = draw_2d(sd_pentagon(uv - (ro + o_pentagon), 0.1), vec3(0.5, 0.5, 0.0), 0.01);
	// draw octogon
	vec2 o_octogon = vec2(1.1, -0.2);
	vec4 c_octogon = draw_2d(sd_octogon(uv - (ro + o_octogon), 0.1), vec3(0.5, 0.5, 0.0), 0.01);
	
	// draw round or annular
	vec2 o_pentagon_round = vec2(0.5, -0.8);
	vec2 o_pentagon_annular = vec2(1.3, -0.8);
	float pentagon_round = sd_round(sd_pentagon(uv - (ro + o_pentagon_round), 0.2), 0.1);
	float pentagon_annular = sd_annular(sd_pentagon(uv - (ro + o_pentagon_annular), 0.2), 0.05);
	vec4 c_pentagon_round = draw_2d(pentagon_round, vec3(0.7, 0.3, 0.2), 0.01);
	vec4 c_pentagon_annular = draw_2d(pentagon_annular, vec3(0.7, 0.3, 0.2), 0.01);

	// draw vesica
	vec2 o_vesica = vec2(1.4, -0.2);
	vec4 c_vesica = draw_2d(sd_vesica(uv - (ro + o_vesica), 0.2, 0.1), vec3(0.5, 0.5, 0.0), 0.01);
	// draw cross
	vec2 o_cross = vec2(1.7, -0.2);
	vec4 c_cross = draw_2d(sd_cross(uv - (ro + o_cross), vec2(0.2, 0.1), 0.05), vec3(0.5, 0.5, 0.0), 0.01);

	// out
	add_color(vec4(1.0, 0.8, 0.7, 1.0));

	add_color(c_line1);
	add_color(c_line2);
	add_color(c_seg_line1);
	add_color(c_seg_line2);
	add_color(c_seg_line3);
	add_color(c_circle);
	add_color(c_circle1);
	add_color(c_circle2);

	add_color(c_circle3);
	add_color(c_ellipse);

	add_color(c_bezier_line);

	add_color(c_rect);
	add_color(c_rhombus);

	add_color(c_triangle1);
	add_color(c_triangle2);
	add_color(c_triangle3);

	add_color(c_trapezoid);
	add_color(c_uneven_capsule);

	add_color(c_pentagon);
	add_color(c_octogon);

	add_color(c_vesica);
	add_color(c_cross);

	add_color(c_pentagon_round);
	add_color(c_pentagon_annular);

	cac_color(frag_color); 

}











