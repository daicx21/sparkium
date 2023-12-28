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
};

#define MATERIAL_TYPE_LAMBERTIAN 0
#define MATERIAL_TYPE_SPECULAR 1
#define MATERIAL_TYPE_TRANSMISSIVE 2
#define MATERIAL_TYPE_PRINCIPLED 3
#define MATERIAL_TYPE_BECKMANN 4
#define MATERIAL_TYPE_EMISSION 5

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

bool SameHemiSphere(vec3 in_direction, vec3 out_direction, vec3 normal_direction) {
  return dot(-in_direction, normal_direction) * dot(out_direction, normal_direction) >= 0;
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
      else {
      }
      return D * cosTheta / (4.0f * dot(-in_direction, wh));
    }
    default: {
      return dot(out_direction, normal_direction) * INV_PI;
    }
  }
}

vec3 bsdf(Material material, vec3 in_direction, vec3 out_direction, vec3 normal_direction) {
  switch (material.material_type) {
    case MATERIAL_TYPE_LAMBERTIAN: {
      if (SameHemiSphere(in_direction, out_direction, normal_direction)) return material.albedo_color * INV_PI;
      return vec3(0.0);
    }
    case MATERIAL_TYPE_PRINCIPLED: {
      float kd = 0.2, ks = 1 - kd;
      int n = 8;
      vec3 spe_out_direction = in_direction - dot(in_direction, normal_direction) * normal_direction * 2.0f;
      float cosine = max(dot(out_direction, spe_out_direction), 0.0f);
      float t1 = kd * INV_PI, t2 = ks * (n+2) * INV_PI / 2.0f * pow(cosine,n);
      return material.albedo_color * (t1+t2);
    }
    case MATERIAL_TYPE_BECKMANN: {
      vec3 wh = normalize((-in_direction + out_direction) * .5);
      float cosTheta = dot(wh, normal_direction), D;
      if (material.alpha == material.beta) {
        float cosTheta2 = cosTheta * cosTheta;
        D = exp(- (1.0f - cosTheta2) / cosTheta2 / (material.alpha * material.alpha))
          * INV_PI / (material.alpha * material.alpha) / (cosTheta2 * cosTheta2);
      }
      else {
      }
      return material.albedo_color * D * Fresnel(material, in_direction, wh) 
        / (1.0f + BeckmannLambda(material, in_direction, wh) + BeckmannLambda(material, out_direction, wh))
        / (4.0f * dot(-in_direction, normal_direction) * dot(out_direction, normal_direction));
    }
    default: {
      return vec3(0.5);
    }
  }
}