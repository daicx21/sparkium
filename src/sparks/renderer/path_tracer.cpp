#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
}

glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) const {
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  HitRecord hit_record;
  const int max_bounce = render_settings_->num_bounces;
  std::mt19937 rd(sample ^ x ^ y);
  std::uniform_real_distribution<> dis(0.0, 1.0);
  for (int i = 0; i < max_bounce; i++) {
    if (dis(rd)>0.9) break;
    throughput/=0.9;
    auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record);
    if (t > 0.0f) {
      auto &material =
          scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
      if (material.material_type == MATERIAL_TYPE_EMISSION) {
        radiance += throughput * material.emission * material.emission_strength;
        break;
      }
      else if (material.material_type == MATERIAL_TYPE_SPECULAR) {
        direction -= 2.0f * dot(hit_record.normal, direction) * hit_record.normal;
        origin = hit_record.position;
      }
      else if (material.material_type == MATERIAL_TYPE_LAMBERTIAN) {
        float a = dis(rd) * PI;
        float b = dis(rd) * PI * 2;
        glm::vec3 hh{sin(b) * cos(a), cos(b) * cos(a), sin(a)};
        if (dot(hh, hit_record.normal) < 0.0f) hh *= -1;
        throughput *= dot(hh, hit_record.normal) * material.albedo_color *
          glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(hit_record.tex_coord)};
        origin = hit_record.position;
        direction = hh;
      }
    } else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(direction)};
      break;
    }
  }
  return radiance;
}
}  // namespace sparks
