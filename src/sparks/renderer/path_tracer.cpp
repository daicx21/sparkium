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
  float area_emission = 130.0f * 105.0f, RR_P = 0.9f;

  glm::vec3 Dir[8];
  float sum[8];

  for (int i = 0;; i++) {
    auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record);
    if (t > 0.0f) {
      auto &material = scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();

      if (dot(hit_record.normal, direction) > 0.0f) hit_record.normal *= -1.0f;

      glm::vec3 hn = hit_record.normal, hp = hit_record.position;
      glm::vec2 ht = hit_record.tex_coord;

      if (material.material_type == MATERIAL_TYPE_EMISSION) {
        radiance += throughput * material.emission * material.emission_strength;
        break;
      }
      else if (material.material_type == MATERIAL_TYPE_SPECULAR) {
        if (dis(rd) > RR_P) break;
        throughput /= RR_P;
        direction -= 2.0f * dot(hit_record.normal, direction) * hit_record.normal;
        origin = hit_record.position;
      }
      else if (material.material_type == MATERIAL_TYPE_LAMBERTIAN) {
        throughput *= (material.albedo_color / PI) *
                      glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(ht)};

        glm::vec3 hhh = glm::vec3(dis(rd) * 130.0f + 213.0f, 548.7f, dis(rd) * 105.0f + 227.0f);
        glm::vec3 dir = normalize(hhh - hp);
        auto tt = scene_->TraceRay(hp, dir, 1.0f, 1e4f, &hit_record);
        auto &material1 = scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
        if (!(dot(dir, hn) <= 0.0f || tt == -1.0f || material1.material_type != MATERIAL_TYPE_EMISSION))
        {
          if (dot(hit_record.normal, dir) > 0.0f) hit_record.normal *= -1.0f;
          radiance += throughput * material1.emission * material1.emission_strength *
                      dot(dir, hn) * (-dot(hit_record.normal, dir)) /
                      dot(hhh - hp, hhh - hp) * area_emission;
        }

        if (dis(rd) > RR_P) break;
        throughput /= RR_P;

        float cnt = 0;

        for (int j = 0; j < 8; j ++)
        {
          float a = acos(sqrt(1 - dis(rd)));
          float b = dis(rd) * PI * 2;
          dir = glm::vec3(sin(b) * cos(a), cos(b) * cos(a), sin(a));
          if (dot(dir, hn) < 0.0f) dir *= -1;
          auto ttt = scene_->TraceRay(hp, dir, 1e-3f, 1e4f, &hit_record);
          auto &material2 = scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
          if (ttt != -1.0f &&
            material2.material_type == MATERIAL_TYPE_EMISSION) { sum[j] = 0; continue; }
          cnt += 1.0f;
          Dir[j] = dir;
          sum[j] = dot(dir, hn);
        }

        if (cnt == 0.0f) break;

        for (int j = 1; j < 8; j ++) sum[j] += sum[j-1];

        float sd = dis(rd) * sum[7];
        int id = 0;
        for (int j = 1; j < 8; j ++) if (sd <= sum[j] && sd > sum[j-1]) id = j;
      
        origin = hp;
        direction = Dir[id];
        throughput *= (PI * 2) * (sum[7] / cnt);
      }
    }
    else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(direction)};
      break;
    }
  }
  return radiance;
}
}  // namespace sparks
