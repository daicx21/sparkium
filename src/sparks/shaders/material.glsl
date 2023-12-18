#include "constants.glsl"

struct Material {
  vec3 albedo_color;
  int albedo_texture_id;
  vec3 emission;
  float emission_strength;
  float refraction_ratio;
  float alpha;
  uint material_type;
};

#define MATERIAL_TYPE_LAMBERTIAN 0
#define MATERIAL_TYPE_SPECULAR 1
#define MATERIAL_TYPE_TRANSMISSIVE 2
#define MATERIAL_TYPE_PRINCIPLED 3
#define MATERIAL_TYPE_EMISSION 4

vec3 bsdf(Material material, vec3 in_direction, vec3 out_direction, vec3 normal_direction) {
  switch (material.material_type) {
    case MATERIAL_TYPE_LAMBERTIAN: {
      if (dot(-in_direction, normal_direction) * dot(out_direction, normal_direction) >= 0) return material.albedo_color * INV_PI;
      return vec3(0.0);
    }
    case MATERIAL_TYPE_PRINCIPLED: {
      float kd = 0.2, ks = 1-kd;
      int n = 8;
      vec3 spe_out_direction = in_direction - dot(in_direction, normal_direction) * normal_direction * 2.0f;
      float cosine = max(dot(out_direction, spe_out_direction), 0.0f);
      float t1 = kd * INV_PI, t2 = ks * (n+2) * INV_PI / 2.0f * pow(cosine,n);
      return material.albedo_color * (t1+t2);
    }
    default: {
      return vec3(0.5);
    }
  }
}