#include "constants.glsl"
#include "random.glsl"

struct Material {
  vec3 albedo_color;
  int albedo_texture_id;
  vec3 emission;
  float emission_strength;
  float refraction_ratio;
  float alpha;
  float beta;
  uint material_type;
  float specular_transmission;
  float metallic;
  float subsurface;
  float specular;
  float roughness;
  float specular_tint;
  float anisotropic;
  float sheen;
  float sheen_tint;
  float clearcoat;
  float clearcoat_gloss;
  float eta;
  int normal_texture_id;
};

#define MATERIAL_TYPE_LAMBERTIAN 0
#define MATERIAL_TYPE_SPECULAR 1
#define MATERIAL_TYPE_TRANSMISSIVE 2
#define MATERIAL_TYPE_PRINCIPLED 3
#define MATERIAL_TYPE_BECKMANN 4
#define MATERIAL_TYPE_EMISSION 5

Material mat;

vec3 win, wout, n, tg1, tg2, ng, h, hl;

float aspect, alpha_x, alpha_y;

float sqr(float x) { return x * x; }

float luminance(vec3 col) {
  return col[0] * 0.212671 + col[1] * 0.715160 + col[2] * 0.072169;
}

vec3 to_local(vec3 w) {
  return vec3(dot(w, tg1), dot(w, tg2), dot(w, n));
}

vec3 to_world(vec3 w) {
  return tg1 * w[0] + tg2 * w[1] + n * w[2];
}

float F_D(vec3 w) {
  return 1.0 + (2.0 * mat.roughness * sqr(abs(dot(h, wout))) - 0.5) * pow(1.0 - abs(dot(n, w)), 5); 
}

float F_SS(vec3 w) {
  return 1.0 + (mat.roughness * sqr(abs(dot(h, wout))) - 1.0) * pow(1.0 - abs(dot(n, w)), 5);
}

float Lambda(vec3 w) {
  vec3 wl = to_local(w);
  return (sqrt(1.0 + (sqr(wl[0] * alpha_x) + sqr(wl[1] * alpha_y)) / (sqr(wl[2]))) - 1.0) * 0.5;
}

float Lambda_c(vec3 w) {
  vec3 wl = to_local(w);
  return (sqrt(1.0 + (sqr(wl[0] * 0.25) + sqr(wl[1] * 0.25)) / (sqr(wl[2]))) - 1.0) * 0.5;
}

float get_Dm() {
  return 1.0 / (PI * alpha_x * alpha_y * sqr(sqr(hl[0] / alpha_x) + sqr(hl[1] / alpha_y) + sqr(hl[2])));
}

float G(vec3 w) {
  return 1.0 / (1.0 + Lambda(w));
}

float get_Dc(float alpha_g) {
  return (sqr(alpha_g) - 1.0) / (PI * 2.0 * log(alpha_g) * (1.0 + (sqr(alpha_g) - 1.0) * sqr(hl[2])));
}

float get_Gc(vec3 w) {
  return 1.0 / (1.0 + Lambda_c(w));
}

bool SameHemiSphere(vec3 in_direction, vec3 out_direction, vec3 normal_direction) {
  return dot(-in_direction, normal_direction) * dot(out_direction, normal_direction) >= 0;
}

vec3 GeneratePerpendicular() {
  vec3 hx;
  if (abs(n.x) <= abs(n.y) && abs(n.x) <= abs(n.z)) {
    hx = vec3(0.0f, -n.z, n.y);
  }
  else if (abs(n.y) <= abs(n.z) && abs(n.y) <= abs(n.x)) {
    hx = vec3(n.z, 0.0f, -n.x);
  }
  else {
    hx = vec3(-n.y, n.x, 0.0f);
  }
  return normalize(hx);
}

mat3 GenerateRotation() {
  vec3 hx = GeneratePerpendicular();
  vec3 hy = cross(n, hx);
  return mat3(hx, hy, n);
}

void pre_gao(Material material, vec3 N, vec3 NG, vec3 TG1) {
  mat = material; n = N; ng = NG;
  tg1 = cross(n, vec3(1.0, 0.0, 0.0));
  if (dot(tg1, tg1) == 0) tg1 = cross(n, vec3(0.0, 1.0, 0.0));
  tg1 = normalize(tg1); tg2 = cross(n, tg1);
}

vec3 sample_diffuse() {
  float a = sqrt(RandomFloat());
  float b = RandomFloat() * PI * 2;
  mat3 rotation = GenerateRotation();
  return rotation * vec3(a * cos(b), a * sin(b), sqrt(1.0f - a * a));
}

vec3 sample_disney(vec3 in_direction) {
  win = -in_direction;
  //return sample_diffuse();
  float alpha_g = (1.0 - mat.clearcoat_gloss) * 0.1 + mat.clearcoat_gloss * 0.001;
  float u_0 = RandomFloat(), u_1 = RandomFloat();
  float phi = acos(0.98);
  float theta = acos(2.0 * PI * u_1);
  h = to_world(vec3(sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi)));
  vec3 res = normalize(2.0 * dot(h, win) * h - win);
  if (dot(res, n) <= 0.0) res = normalize(2.0 * dot(n, win) * n - win);
  return res;
}

vec3 sample_bsdf(vec3 in_direction)
{
  switch (mat.material_type) {
    case MATERIAL_TYPE_PRINCIPLED: {
      return sample_disney(in_direction);
    }
    default: {
      return sample_diffuse();
    }
  }
}

float pdf_diffuse(vec3 out_direction) {
  return max(dot(out_direction, n), 0) * INV_PI;
}

float pdf_disney(vec3 in_direction, vec3 out_direction) {
  win = -in_direction; wout = out_direction;
  h = normalize(win + wout); hl = to_local(h);
  //return pdf_diffuse(out_direction);
  float alpha_g = (1.0 - mat.clearcoat_gloss) * 0.1 + mat.clearcoat_gloss * 0.001;
  float Dc = get_Dc(alpha_g);
  return Dc * abs(dot(n, h)) / (4.0 * abs(dot(h, wout)));
}

float pdf(vec3 in_direction, vec3 out_direction) {
  switch (mat.material_type) {
    case MATERIAL_TYPE_PRINCIPLED: {
      return pdf_disney(in_direction, out_direction);
    }
    default: {
      return pdf_diffuse(out_direction);
    }
  }
}

vec3 bsdf_disney(vec3 in_direction, vec3 out_direction) {
  win = -in_direction; wout = out_direction;
  h = normalize(win + wout); hl = to_local(h);
  
  vec3 f_baseDiffuse = mat.albedo_color * INV_PI * F_D(win) * F_D(wout) * abs(dot(n, wout));
  vec3 f_subsurface = mat.albedo_color * 1.25 * INV_PI *
                      (F_SS(win) * F_SS(wout) * (1.0 / (abs(dot(n, win)) + abs(dot(n, wout))) - 0.5) + 0.5) * abs(dot(n, wout));
  vec3 f_diffuse = (1.0 - mat.subsurface) * f_baseDiffuse + mat.subsurface * f_subsurface;
  
  vec3 Fm = mat.albedo_color + (vec3(1.0) - mat.albedo_color) * pow(1.0 - abs(dot(h, wout)), 5);
  aspect = sqrt(1.0 - 0.9 * mat.anisotropic);
  alpha_x = max(0.0001, sqr(mat.roughness) / aspect);
  alpha_y = max(0.0001, sqr(mat.roughness) * aspect);
  float Dm = get_Dm();
  float Gm = G(win) * G(wout);
  vec3 f_metal = Fm * Dm * Gm / (4.0 * abs(dot(n, win)));

  vec3 Fc = vec3(0.04 + 0.96 * pow(1.0 - abs(dot(h, wout)), 5));
  float alpha_g = (1.0 - mat.clearcoat_gloss) * 0.1 + mat.clearcoat_gloss * 0.001;
  float Dc = get_Dc(alpha_g);
  float Gc = get_Gc(win) * get_Gc(wout);
  vec3 f_clearcoat = Fc * Dc * Gc / (4.0 * abs(dot(n, win)));

  /*
  float lum = luminance(mat.albedo_color);
  vec3 C_tint = vec3(1.0);
  if (lum > 0) C_tint = mat.albedo_color / lum;
  vec3 C_sheen = vec3(1.0 - mat.sheen_tint) + mat.sheen_tint * C_tint;
  vec3 f_sheen = C_sheen * pow(1.0 - abs(dot(h, wout)), 5) * abs(dot(n, wout));*/

  return f_clearcoat;
}

vec3 bsdf(vec3 in_direction, vec3 out_direction) {
  switch (mat.material_type) {
    case MATERIAL_TYPE_LAMBERTIAN: {
      if (SameHemiSphere(in_direction, out_direction, n))
        return mat.albedo_color * dot(out_direction, n) * INV_PI;
      return vec3(0.0);
    }
    case MATERIAL_TYPE_PRINCIPLED: {
      return min(bsdf_disney(in_direction, out_direction), vec3(1.0f));
    }
    default: {
      return vec3(0.5);
    }
  }
}