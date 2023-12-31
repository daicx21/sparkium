#include "constants.glsl"
#include "random.glsl"

struct Material {
  vec3 albedo_color;
  int albedo_texture_id;
  vec3 emission;
  float emission_strength;
  int normal_texture_id;
  float alpha;
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
};

#define MATERIAL_TYPE_LAMBERTIAN 0
#define MATERIAL_TYPE_SPECULAR 1
#define MATERIAL_TYPE_TRANSMISSIVE 2
#define MATERIAL_TYPE_PRINCIPLED 3
#define MATERIAL_TYPE_EMISSION 4

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

vec3 gao_1, gao_2;

void gao_frame(vec3 n) {
  if (n[2] < -1.0 + eps) {
    gao_1 = vec3(0.0, -1.0, 0.0);
    gao_2 = vec3(-1.0, 0.0, 0.0);
  }
  else {
    float a = 1.0 / (1.0 + n[2]), b = -n[0] * n[1] * a;
    gao_1 = vec3(1.0 - n[0] * n[0] * a, b, -n[0]);
    gao_2 = vec3(b, 1.0 - n[1] * n[1] * a, -n[1]);
  }
}

float R0(float eta) {
  return sqr(eta - 1) / sqr(eta + 1);
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
  float a2 = sqr(alpha_g);
  return (a2 - 1.0) / (PI * log(max(a2, eps)) * (1.0 + (a2 - 1.0) * sqr(hl[2])));
}

float get_Gc(vec3 w) {
  return 1.0 / (1.0 + Lambda_c(w));
}

float fresnel_dielectric1(float n_dot_i, float n_dot_t, float eta) {
    float rs = (n_dot_i - eta * n_dot_t) / (n_dot_i + eta * n_dot_t);
    float rp = (eta * n_dot_i - n_dot_t) / (eta * n_dot_i + n_dot_t);
    return (rs * rs + rp * rp) * 0.5;
}

float fresnel_dielectric(float n_dot_i, float eta) {
  float n_dot_t_sq = 1.0 - (1.0 - n_dot_i * n_dot_i) / (eta * eta);
  if (n_dot_t_sq < 0.0) return 1.0;
  float n_dot_t = sqrt(n_dot_t_sq);
  return fresnel_dielectric1(abs(n_dot_i), n_dot_t, eta);
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
  tg1 = TG1; tg2 = cross(n, tg1);
}

vec3 sample_visible_normals_aniso(vec3 wl) {
  float flag = 1.0;
  if (wl[2] < 0.0) wl = -wl, flag = -1.0;
  vec3 hemi_dir_in = normalize(vec3(alpha_x * wl[0], alpha_y * wl[1], wl[2]));
  float u_0 = RandomFloat(), u_1 = RandomFloat();
  float r = sqrt(u_0), phi = 2.0 * PI * u_1;
  float t1 = r * cos(phi), t2 = r * sin(phi);
  float s = (1.0 + hemi_dir_in.z) * 0.5;
  t2 = (1.0 - s) * sqrt(1.0 - t1 * t1) + s * t2;
  gao_frame(hemi_dir_in);
  vec3 hemi_N = t1 * gao_1 + t2 * gao_2 + sqrt(max(0.0, 1.0 - t1 * t1 - t2 * t2)) * hemi_dir_in;
  return normalize(vec3(alpha_x * hemi_N.x, alpha_y * hemi_N.y, max(0.0, hemi_N.z))) * flag;
}

vec3 sample_diffuse() {
  float a = sqrt(RandomFloat());
  float b = RandomFloat() * PI * 2;
  mat3 rotation = GenerateRotation();
  return rotation * vec3(a * cos(b), a * sin(b), sqrt(1.0f - a * a));
}

vec3 sample_metal() {
  aspect = sqrt(1.0 - 0.9 * mat.anisotropic);
  alpha_x = max(0.0001, sqr(mat.roughness) / aspect);
  alpha_y = max(0.0001, sqr(mat.roughness) * aspect);
  h = to_world(sample_visible_normals_aniso(to_local(win)));
  return normalize(2.0 * dot(win, h) * h - win);
}

vec3 sample_clearcoat() {
  float alpha_g = (1.0 - mat.clearcoat_gloss) * 0.1 + mat.clearcoat_gloss * 0.004, a2 = sqr(alpha_g);
  float u_0 = RandomFloat(), u_1 = RandomFloat();
  float cos_phi = sqrt((1.0 - pow(a2, 1.0 - u_0)) / (1.0 - a2)), sin_phi = sqrt(1.0 - sqr(cos_phi));
  float theta = 2.0 * PI * u_1;
  h = to_world(vec3(sin_phi * cos(theta), sin_phi * sin(theta), cos_phi));
  return 2.0 * dot(h, win) * h - win;
}

vec3 sample_glass() {
  float eta = ((dot(ng, win) > 0.0) ? (mat.eta) : (1.0 / mat.eta));
  aspect = sqrt(1.0 - 0.9 * mat.anisotropic);
  alpha_x = max(0.0001, sqr(mat.roughness) / aspect);
  alpha_y = max(0.0001, sqr(mat.roughness) * aspect);
  h = to_world(sample_visible_normals_aniso(to_local(win)));
  if (dot(h, n) < 0.0) h = -h;
  float h_dot_in = dot(h, win);
  float F = fresnel_dielectric(h_dot_in, eta);
  if (fract(RandomFloat()) <= F) return normalize(2.0 * dot(win, h) * h - win);
  float h_dot_out_sq = 1.0 - (1.0 - h_dot_in * h_dot_in) / (eta * eta);
  if (h_dot_out_sq < 0.0) h_dot_out_sq = 0.0;
  if (h_dot_in < 0) h = -h;
  float h_dot_out = sqrt(h_dot_out_sq);
  return normalize(-win / eta + (abs(h_dot_in) / eta - h_dot_out) * h);
}

vec3 sample_disney() {
  if (dot(n, win) * dot(ng, win) < 0.0) n = -n, tg2 = -tg2;
  if (dot(n, win) <= 0.0) return sample_glass();
  float sw_diffuse = (1.0 - mat.specular_transmission) * (1.0 - mat.metallic);
  float sw_metal = (1.0 - mat.specular_transmission * (1.0 - mat.metallic));
  float sw_clearcoat = 0.25 * mat.clearcoat;
  float sw_glass = (1.0 - mat.metallic) * mat.specular_transmission;
  float cdf_clearcoat = sw_clearcoat;
  float cdf_metal = cdf_clearcoat + sw_metal;
  float cdf_glass = cdf_metal+ sw_glass;
  float u_0 = fract(RandomFloat());
  if (u_0 < cdf_clearcoat) return sample_clearcoat();
  if (u_0 < cdf_metal) return sample_metal();
  if (u_0 < cdf_glass) return sample_glass();
  return sample_diffuse();
}

vec3 sample_bsdf(vec3 in_direction)
{
  win = -in_direction;
  switch (mat.material_type) {
    case MATERIAL_TYPE_PRINCIPLED: {
      return sample_disney();
    }
    default: {
      return sample_diffuse();
    }
  }
}

float pdf_diffuse() {
  if (dot(win, ng) <= 0 || dot(wout, ng) <= 0) return 0;
  return max(dot(wout, n), 0) * INV_PI;
}

float pdf_metal() {
  if (dot(win, ng) <= 0 || dot(wout, ng) <= 0) return 0;
  if (dot(n, wout) <= 0 || dot(n, h) <= 0) return 0;
  aspect = sqrt(1.0 - 0.9 * mat.anisotropic);
  alpha_x = max(0.0001, sqr(mat.roughness) / aspect);
  alpha_y = max(0.0001, sqr(mat.roughness) * aspect);
  float Dm = get_Dm();
  float Gm = G(win) * G(wout);
  return (Gm * Dm) / (4.0 * dot(n, win));
}

float pdf_clearcoat() {
  if (dot(win, ng) <= 0.0f || dot(wout, ng) <= 0.0f) return 0.0;
  float alpha_g = (1.0 - mat.clearcoat_gloss) * 0.1 + mat.clearcoat_gloss * 0.004;
  return get_Dc(alpha_g) * abs(dot(n, h)) / (4.0 * abs(dot(h, wout)));
}

float pdf_glass() {
  bool reflect = (dot(ng, win) * dot(ng, wout) > 0);
  float eta = ((dot(ng, win) > 0) ? (mat.eta) : (1.0 / mat.eta));
  if (!reflect) h = -normalize(win + wout * eta), hl = to_local(h);
  aspect = sqrt(1.0 - 0.9 * mat.anisotropic);
  alpha_x = max(0.0001, sqr(mat.roughness) / aspect);
  alpha_y = max(0.0001, sqr(mat.roughness) * aspect);
  float h_dot_in = dot(h, win);
  float Fg = fresnel_dielectric(dot(h, win), eta), Dg = get_Dm(), Gg = G(win) * G(wout);
  if (reflect) return (Fg * Dg * Gg) / (4.0 * abs(dot(n, win)));
  float h_dot_out = dot(h, wout);
  float sqrt_denom = h_dot_in + eta * h_dot_out;
  float dh_dout = eta * eta * h_dot_out / sqr(sqrt_denom);
  return (1.0 - Fg) * Dg * Gg * abs(dh_dout * h_dot_in / dot(n, win));
}

float pdf_disney() {
  if (dot(n, win) * dot(ng, win) < 0) n = -n, tg2 = -tg2, hl = to_local(h);
  float sw_diffuse = (1.0 - mat.specular_transmission) * (1.0 - mat.metallic);
  float sw_metal = (1.0 - mat.specular_transmission * (1.0 - mat.metallic));
  float sw_clearcoat = 0.25 * mat.clearcoat;
  float sw_glass = (1.0 - mat.metallic) * mat.specular_transmission;
  if (dot(n, win) < 0 || dot(n, wout) < 0) return pdf_glass();
  return sw_clearcoat * pdf_clearcoat() + sw_metal * pdf_metal() + sw_diffuse * pdf_diffuse() + sw_glass * pdf_glass();
}

float pdf(vec3 in_direction, vec3 out_direction) {
  win = -in_direction; wout = out_direction;
  h = normalize(win + wout); hl = to_local(h);
  switch (mat.material_type) {
    case MATERIAL_TYPE_PRINCIPLED: {
      return pdf_disney();
    }
    default: {
      return pdf_diffuse();
    }
  }
}

vec3 bsdf_disney(vec3 in_direction, vec3 out_direction) {
  win = -in_direction; wout = out_direction;
  h = normalize(win + wout); hl = to_local(h);
  if (dot(n, win) * dot(ng, win) < 0) n = -n, tg2 = -tg2, hl = to_local(h);
  
  vec3 f_baseDiffuse = mat.albedo_color * INV_PI * F_D(win) * F_D(wout) * abs(dot(n, wout));
  vec3 f_subsurface = mat.albedo_color * 1.25 * INV_PI *
                      (F_SS(win) * F_SS(wout) * (1.0 / (abs(dot(n, win)) + abs(dot(n, wout))) - 0.5) + 0.5) * abs(dot(n, wout));
  vec3 f_diffuse = (1.0 - mat.subsurface) * f_baseDiffuse + mat.subsurface * f_subsurface;
  
  float eta = ((dot(ng, win) > 0) ? (mat.eta) : (1.0 / mat.eta));
  float lum = luminance(mat.albedo_color);
  vec3 C_tint = vec3(1.0);
  if (lum > 0) C_tint = mat.albedo_color / lum;
  vec3 Ks = vec3(1.0 - mat.specular_tint) + mat.specular_tint * C_tint;
  vec3 C0 = mat.specular * (R0(eta)) * (1 - mat.metallic) * Ks + mat.metallic * mat.albedo_color;
  vec3 Fm = C0 + (vec3(1.0) - C0) * pow(1.0 - abs(dot(h, wout)), 5);
  aspect = sqrt(1.0 - 0.9 * mat.anisotropic);
  alpha_x = max(0.0001, sqr(mat.roughness) / aspect);
  alpha_y = max(0.0001, sqr(mat.roughness) * aspect);
  float Dm = get_Dm();
  float Gm = G(win) * G(wout);
  vec3 f_metal = Fm * Dm * Gm / (4.0 * abs(dot(n, win)));

  vec3 Fc = vec3(0.04 + 0.96 * pow(1.0 - abs(dot(h, wout)), 5));
  float alpha_g = (1.0 - mat.clearcoat_gloss) * 0.1 + mat.clearcoat_gloss * 0.004;
  float Dc = get_Dc(alpha_g);
  float Gc = get_Gc(win) * get_Gc(wout);
  vec3 f_clearcoat = Fc * Dc * Gc / (4.0 * abs(dot(n, win)));

  bool reflect = (dot(ng, win) * dot(ng, wout) > 0);
  if (!reflect) h = -normalize(win + wout * eta), hl = to_local(h);
  float Fg = fresnel_dielectric(dot(h, win), eta), Dg = get_Dm(), Gg = G(win) * G(wout);
  vec3 f_glass;
  if (reflect) f_glass = mat.albedo_color * Fg * Dg * Gg / (4.0 * abs(dot(n, win)));
  else f_glass = sqrt(mat.albedo_color) * (1 - Fg) * Dg * Gg * abs(dot(h, win) * dot(h, wout)) /
                 (abs(dot(n, win)) * sqr(dot(h, win) + eta * dot(h, wout)));

  vec3 C_sheen = vec3(1.0 - mat.sheen_tint) + mat.sheen_tint * C_tint;
  vec3 f_sheen = C_sheen * pow(1.0 - abs(dot(h, wout)), 5) * abs(dot(n, wout));

  float w_diffuse = (1.0 - mat.specular_transmission) * (1.0 - mat.metallic);
  float w_sheen = (1.0 - mat.metallic) * mat.sheen;
  float w_metal = (1.0 - mat.specular_transmission * (1.0 - mat.metallic));
  float w_clearcoat = 0.25 * mat.clearcoat;
  float w_glass = (1.0 - mat.metallic) * mat.specular_transmission;

  vec3 f_disney;
  if (dot(n, win) < 0 || dot(n, wout) < 0) f_disney = w_glass * f_glass;
  else
  {
    if (reflect) f_disney = w_diffuse * f_diffuse + w_sheen * f_sheen + w_metal * f_metal + w_clearcoat * f_clearcoat + w_glass * f_glass;
    else f_disney = w_glass * f_glass;
  }
  return f_disney;
}

vec3 bsdf(vec3 in_direction, vec3 out_direction) {
  switch (mat.material_type) {
    case MATERIAL_TYPE_LAMBERTIAN: {
      if (SameHemiSphere(in_direction, out_direction, n))
        return mat.albedo_color * dot(out_direction, n) * INV_PI;
      return vec3(0.0);
    }
    case MATERIAL_TYPE_PRINCIPLED: {
      return bsdf_disney(in_direction, out_direction);
    }
    default: {
      return vec3(0.5);
    }
  }
}