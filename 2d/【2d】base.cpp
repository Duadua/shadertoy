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

float sd_rect(vec2 p, vec2 size) {
	vec2 d = abs(p) - size;										// abs(vec2) -- 每个分量取绝对值
	return length(max(d, vec2(0.0))) + min(max(d.x, d.y), 0.0);	// max(vec2, vec2) -- 每个分量取最大值 
	// 取差向量的 正分量向量长度 + 最大的负分量
}

float sd_rhombus(vec2 p, vec2 b) {							// 菱形
	vec2 q = abs(p);
	
	vec2 c = 0.5*b;
	vec2 e = q - c;
	c.x = -1.0*c.x;
	float h = clamp(dot(e, c) / dot(c, c), -1.0, 1.0);
	float d = length(e - c*h);
	return d*sign(cross2(e, c));
}

// ============================================================

vec4 draw_2d(float sd, vec3 color, float att) {
	float t = smoothstep(0.0, att, sd);
	return vec4(color, 1.0 - t);
}

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

	// draw point (circle)
	vec2 o = vec2(0.0, 0.0);
	vec2 o1 = vec2(0.0, 1.0);
	vec2 o2 = vec2(1.0, 0.0);
	vec4 c_circle = draw_2d(sd_circle(uv - (ro + o), 0.02), vec3(0.4,0.2,0.6), 0.01);
	vec4 c_circle1 = draw_2d(sd_circle(uv - (ro + o1), 0.02), vec3(0.4,0.2,0.6), 0.01);
	vec4 c_circle2 = draw_2d(sd_circle(uv - (ro + o2), 0.02), vec3(0.4,0.2,0.6), 0.01);

	// draw line
	vec4 c_line1 = draw_2d(sd_line(uv, ro + o, ro + o1, 0.01), vec3(0.5), 0.01);
	vec4 c_line2 = draw_2d(sd_line(uv, ro + o, ro + o2, 0.01), vec3(0.5), 0.01);
	// draw segment line
	vec4 c_seg_line1 = draw_2d(sd_seg_line(uv, ro + o, ro + o1, 0.01), vec3(0.5, 1.0, 0.2), 0.01);
	vec4 c_seg_line2 = draw_2d(sd_seg_line(uv, ro + o, ro + o2, 0.01), vec3(0.5, 1.0, 0.2), 0.01);
	vec4 c_seg_line3 = draw_2d(sd_seg_line(uv, ro + o1, ro + o2, 0.01), vec3(0.5, 1.0, 0.2), 0.01);
	
	// draw_triangle
	vec2 B = vec2(-0.9, 0.8);
	vec2 A = vec2(-0.5, 0.5);
	vec2 C = vec2(-0.7, -0.3);
	vec4 c_triangle = draw_triangle(uv, ro + o, A, B, C, vec4(0.0, 1.0, 0.0, 1.0));


	// draw rect 
	vec2 o_rect = vec2(0.6, -0.4);
	vec4 c_rect = draw_2d(sd_rect(uv - (ro + o_rect), vec2(0.5, 0.3)), vec3(1.0, 0.0, 0.0), 0.01);

	// draw rhombus
	vec2 o_rhombus = vec2(-0.9, 0.9);
	vec4 c_rhombus = draw_2d(sd_rhombus(uv - (ro + o_rhombus), vec2(0.1, 0.2)), vec3(0.0, 0.0, 1.0), 0.01);

	// out
	add_color(vec4(1.0, 0.8, 0.7, 1.0));

	//add_color(c_triangle);
	add_color(c_rect);
	add_color(c_rhombus);

	add_color(c_line1);
	add_color(c_line2);
	add_color(c_seg_line1);
	add_color(c_seg_line2);
	add_color(c_seg_line3);
	add_color(c_circle);
	add_color(c_circle1);
	add_color(c_circle2);

	cac_color(frag_color); 

}




