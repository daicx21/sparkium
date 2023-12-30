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
};

#define MATERIAL_TYPE_LAMBERTIAN 0
#define MATERIAL_TYPE_SPECULAR 1
#define MATERIAL_TYPE_TRANSMISSIVE 2
#define MATERIAL_TYPE_PRINCIPLED 3
#define MATERIAL_TYPE_BECKMANN 4
#define MATERIAL_TYPE_EMISSION 5

float sqr(float x) { return x * x; }

float Fresnel(Material material, vec3 in_direction, vec3 normal_direction) {
  float eta = material.refraction_ratio;
  float cos_t = dot(in_direction, normal_direction), sin_t = sqrt(1.0f - cos_t * cos_t);
  if (cos_t < 0.0f) cos_t = -cos_t;
  else eta = 1.0f / eta;
  float sin_p = sin_t / eta, cos_p, fresnel;
  if (sin_p > 1) fresnel = 1.0f;
  else
  {
    cos_p = sqrt(1.0f - sin_p * sin_p);
    fresnel = (pow(((cos_t - eta * cos_p) / (cos_t + eta * cos_p)), 2.0f) +
      pow(((cos_p - eta * cos_t) / (cos_p + eta * cos_t)), 2.0f)) / 2.0f;
  }
  return fresnel;
}

float BeckmannLambda(Material material, vec3 direction, vec3 normal) {
  float cosTheta = dot(direction, normal);
  float tanTheta = abs(sqrt(max(0.0f, 1.0f - cosTheta * cosTheta)) / cosTheta);
  float a = 1.0f / (material.alpha * tanTheta);
  if (a >= 1.6f) return 0;
  else return (1.0f - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
}

bool SameHemiSphere(vec3 in_direction, vec3 out_direction, vec3 normal_direction) {
  return dot(-in_direction, normal_direction) * dot(out_direction, normal_direction) >= 0;
}

vec3 GeneratePerpendicular(vec3 normal) {
  vec3 hx;
  if (abs(normal.x) <= abs(normal.y) && abs(normal.x) <= abs(normal.z)) {
    hx = vec3(0.0f, -normal.z, normal.y);
  }
  else if (abs(normal.y) <= abs(normal.z) && abs(normal.y) <= abs(normal.x)) {
    hx = vec3(normal.z, 0.0f, -normal.x);
  }
  else {
    hx = vec3(-normal.y, normal.x, 0.0f);
  }
  return normalize(hx);
}

mat3 GenerateRotation(vec3 normal) {
  vec3 hx = GeneratePerpendicular(normal);
  vec3 hy = cross(normal, hx);
  return mat3(hx, hy, normal);
}

vec3 SampleWh(Material material, mat3 rotation, vec3 in_direction) {
  if (material.alpha == material.beta) {
    float phi = 2 * PI * RandomFloat();
    float logFloat = log(RandomFloat());
    if (isinf(logFloat)) logFloat = 0;
    float tanTheta2 = -material.alpha * material.alpha * logFloat;
    float cosTheta = 1.0f / sqrt(1.0f + tanTheta2), sinTheta = sqrt(max(0.0f, 1.0f - cosTheta * cosTheta));
    vec3 wh = rotation * vec3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
    if (!SameHemiSphere(in_direction, wh, rotation[2])) wh *= -1;
    return wh;
  }
  else {
    /*float u = RandomFloat();
    float phi = atan(material.beta / material.alpha * tan(2 * PI * u + 0.5f * PI));
    if (u > 0.5f) phi += PI;
    u = log(RandomFloat());
    float tanTheta2 = -u / (cos(phi) * cos(phi) / (material.alpha * material.alpha) 
      + sin(phi) * sin(phi) / (material.beta * material.beta));
    float cosTheta = 1.0f / sqrt(1.0f + tanTheta2), sinTheta = sqrt(max(0.0f, 1.0f - cosTheta * cosTheta));
    vec3 wh = rotation * vec3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
    if (!SameHemiSphere(in_direction, wh, rotation[2])) wh *= -1;
    return wh;*/
  }
}

vec3 sample_beckmann(Material material, vec3 in_direction, vec3 normal_direction) {
  mat3 rotation = GenerateRotation(normal_direction);
  vec3 wh = SampleWh(material, rotation, in_direction);
  return in_direction - 2.0f * dot(wh, in_direction) * wh;
}

vec3 sample_diffuse(vec3 normal_direction) {
  float a = sqrt(RandomFloat());
  float b = RandomFloat() * PI * 2;
  mat3 rotation = GenerateRotation(normal_direction);
  return rotation * vec3(a * cos(b), a * sin(b), sqrt(1.0f - a * a));
}

vec3 sample_disney(Material material, vec3 in_direction, vec3 normal_direction) {
  return sample_diffuse(normal_direction);
}

vec3 sample_bsdf(Material material, vec3 in_direction, vec3 normal_direction)
{
  switch (material.material_type) {
    case MATERIAL_TYPE_BECKMANN: {
      return sample_beckmann(material, in_direction, normal_direction);
    }
    case MATERIAL_TYPE_PRINCIPLED: {
      return sample_disney(material, in_direction, normal_direction);
    }
    default: {
      return sample_diffuse(normal_direction);
    }
  }
}

float pdf_diffuse(vec3 out_direction, vec3 normal_direction) {
  return max(dot(out_direction, normal_direction), 0) * INV_PI;
}

float pdf_disney(vec3 in_direction, vec3 out_direction, vec3 normal_direction) {
  return pdf_diffuse(out_direction, normal_direction);
}

float pdf(Material material, vec3 in_direction, vec3 out_direction, vec3 normal_direction) {
  switch (material.material_type) {
    case MATERIAL_TYPE_BECKMANN: {
      vec3 wh = normalize((-in_direction + out_direction) * .5);
      float cosTheta = dot(wh, normal_direction), D;
      if (material.alpha == material.beta) {
        float cosTheta2 = cosTheta * cosTheta;
        D = exp(- (1.0f - cosTheta2) / cosTheta2 / (material.alpha * material.alpha))
          * INV_PI / (material.alpha * material.alpha) / (cosTheta2 * cosTheta2);
      }
      return D * cosTheta / (4.0f * dot(-in_direction, wh));
    }
    case MATERIAL_TYPE_PRINCIPLED: {
      return pdf_disney(in_direction, out_direction, normal_direction);
    }
    default: {
      return pdf_diffuse(out_direction, normal_direction);
    }
  }
}

Material mat;

vec3 win, wout, n, ng, h, hl;

float aspect, alpha_x, alpha_y;

float luminance(vec3 col) {
  return col[0] * 0.212671 + col[1] * 0.715160 + col[2] * 0.072169;
}

vec3 to_local(vec3 w) {
  float z = dot(w, n), hh = 1.0 - sqr(z);
  float x = sqrt(hh * 0.5);
  return vec3(x,x,z);
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

float G(vec3 w) {
  return 1.0 / (1.0 + Lambda(w));
}

float G_c(vec3 w) {
  return 1.0 / (1.0 + Lambda_c(w));
}

vec3 bsdf_disney(Material material, vec3 in_direction, vec3 out_direction, vec3 normal_direction) {
  mat = material; win = -in_direction; wout = out_direction; n = normal_direction; ng = n;
  h = normalize(win + wout); hl = to_local(h);
  
  vec3 f_baseDiffuse = mat.albedo_color * INV_PI * F_D(in_direction) * F_D(out_direction) * abs(dot(n, wout));
  vec3 f_subsurface = mat.albedo_color * 1.25 * INV_PI *
                      (F_SS(win) * F_SS(wout) * (1.0 / (abs(dot(n, win)) + abs(dot(n, wout))) - 0.5) + 0.5) * abs(dot(n, wout));
  vec3 f_diffuse = (1.0 - mat.subsurface) * f_baseDiffuse + mat.subsurface * f_subsurface;
  
  vec3 F_m = mat.albedo_color + (vec3(1.0) - mat.albedo_color) * pow(1.0 - abs(dot(h, wout)), 5);
  aspect = sqrt(1.0 - 0.9 * mat.anisotropic);
  alpha_x = max(0.0001, sqr(mat.roughness) / aspect);
  alpha_y = max(0.0001, sqr(mat.roughness) * aspect);
  float D_m = 1.0 / (PI * alpha_x * alpha_y * sqr(sqr(hl[0] / alpha_x) + sqr(hl[1] / alpha_y) + sqr(hl[2])));
  float G_m = G(win) * G(wout);
  vec3 f_metal = F_m * D_m * G_m / (4.0 * abs(dot(n, win)));

  vec3 F_c = vec3(0.04 + 0.96 * pow(1.0 - abs(dot(h, wout)), 5));
  float alpha_g = (1.0 - mat.clearcoat_gloss) * 0.1 + mat.clearcoat_gloss * 0.001;
  float D_c = (sqr(alpha_g) - 1.0) / (PI * 2.0 * log(alpha_g) * (1.0 + (sqr(alpha_g) - 1.0) * sqr(hl[2])));
  float G_c = G_c(win) * G_c(wout);
  vec3 f_clearcoat = F_c * D_c * G_c / (4.0 * abs(dot(n, win)));

  vec3 f_glass;

  if (dot(ng, win) * dot(ng, wout) > 0) f_glass = mat.albedo_color * F_g * D_g * G_g / (4.0 * abs(dot(n, win)));

  float lum = luminance(mat.albedo_color);
  vec3 C_tint = vec3(1.0);
  if (lum > 0) C_tint = mat.albedo_color / lum;
  vec3 C_sheen = vec3(1.0 - mat.sheen_tint) + mat.sheen_tint * C_tint;
  vec3 f_sheen = C_sheen * pow(1.0 - abs(dot(h, wout)), 5) * abs(dot(n, wout));

  return f_metal;
}

vec3 bsdf(Material material, vec3 in_direction, vec3 out_direction, vec3 normal_direction) {
  switch (material.material_type) {
    case MATERIAL_TYPE_LAMBERTIAN: {
      if (SameHemiSphere(in_direction, out_direction, normal_direction))
        return material.albedo_color * dot(out_direction, normal_direction) * INV_PI;
      return vec3(0.0);
    }
    case MATERIAL_TYPE_PRINCIPLED: {
      return min(bsdf_disney(material, in_direction, out_direction, normal_direction), vec3(1.0f));
    }
    case MATERIAL_TYPE_BECKMANN: {
      vec3 wh = normalize((-in_direction + out_direction) * .5);
      float cosTheta = dot(wh, normal_direction), D;
      if (material.alpha == material.beta) {
        float cosTheta2 = cosTheta * cosTheta;
        D = exp(- (1.0f - cosTheta2) / cosTheta2 / (material.alpha * material.alpha))
          * INV_PI / (material.alpha * material.alpha) / (cosTheta2 * cosTheta2);
      }
      return material.albedo_color * D * Fresnel(material, in_direction, wh) 
        / (1.0f + BeckmannLambda(material, in_direction, wh) + BeckmannLambda(material, out_direction, wh))
        / (4.0f * dot(-in_direction, normal_direction) * dot(out_direction, normal_direction)) * dot(out_direction, normal_direction);
    }
    default: {
      return vec3(0.5);
    }
  }
}