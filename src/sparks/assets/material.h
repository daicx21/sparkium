#pragma once
#include "cstdint"
#include "glm/glm.hpp"
#include "sparks/assets/util.h"

namespace sparks {

enum MaterialType : int {
  MATERIAL_TYPE_LAMBERTIAN = 0,
  MATERIAL_TYPE_SPECULAR = 1,
  MATERIAL_TYPE_TRANSMISSIVE = 2,
  MATERIAL_TYPE_PRINCIPLED = 3,
  MATERIAL_TYPE_EMISSION = 4
};

class Scene;

struct Material {
  glm::vec3 albedo_color{0.8f};
  int albedo_texture_id{0};
  glm::vec3 emission{0.0f};
  float emission_strength{1.0f};
  int normal_texture_id{-1};
  float alpha{1.0f};
  MaterialType material_type{MATERIAL_TYPE_LAMBERTIAN};
  float specular_transmission{0.5f};
  float metallic{0.5f};
  float subsurface{0.5f};
  float specular{0.5f};
  float roughness{0.5f};
  float specular_tint{0.5f};
  float anisotropic{0.5f};
  float sheen{0.5f};
  float sheen_tint{0.5f};
  float clearcoat{0.5f};
  float clearcoat_gloss{0.5f};
  float eta{1.5f};

  float reserve[1]{};
  Material() = default;
  explicit Material(const glm::vec3 &albedo);
  Material(Scene *scene, const tinyxml2::XMLElement *material_element);
};
}  // namespace sparks
