#include "constants.glsl"

struct HitRecord {
  int hit_entity_id;
  vec3 position;
  vec3 normal;
  vec3 geometry_normal;
  vec3 tangent;
  vec2 tex_coord;
};

int inside_id;

vec3 gao1(vec3 n) {
  if (n[2] < -1.0 + eps) return vec3(0.0, -1.0, 0.0);
  float a = 1.0 / (1.0 + n[2]), b = -n[0] * n[1] * a;
  return vec3(1.0 - n[0] * n[0] * a, b, -n[0]);
}

HitRecord GetHitRecord(RayPayload ray_payload, vec3 origin, vec3 direction) {
  HitRecord hit_record;
  ObjectInfo object_info = object_infos[ray_payload.object_id];
  Vertex v0 = GetVertex(
      object_info.vertex_offset +
      indices[object_info.index_offset + ray_payload.primitive_id * 3 + 0]);
  Vertex v1 = GetVertex(
      object_info.vertex_offset +
      indices[object_info.index_offset + ray_payload.primitive_id * 3 + 1]);
  Vertex v2 = GetVertex(
      object_info.vertex_offset +
      indices[object_info.index_offset + ray_payload.primitive_id * 3 + 2]);
  hit_record.hit_entity_id = int(ray_payload.object_id);

  mat3 object_to_world = mat3(ray_payload.object_to_world);
  hit_record.position = ray_payload.object_to_world *
                        vec4(mat3(v0.position, v1.position, v2.position) *
                                 ray_payload.barycentric,
                             1.0);

  hit_record.normal = normalize(transpose(inverse(object_to_world)) *
                                mat3(v0.normal, v1.normal, v2.normal) *
                                ray_payload.barycentric);
  hit_record.geometry_normal =
      normalize(transpose(inverse(object_to_world)) *
                cross(v1.position - v0.position, v2.position - v0.position));
  hit_record.tangent = object_to_world * mat3(v0.tangent, v1.tangent, v2.tangent) * ray_payload.barycentric;
  hit_record.tangent = normalize(hit_record.tangent
      - dot(hit_record.tangent, hit_record.normal) * hit_record.normal);

  if (dot(direction, hit_record.geometry_normal) > 0.0) hit_record.geometry_normal = -hit_record.geometry_normal;
  if (hit_record.hit_entity_id == inside_id && inside_id != 0) hit_record.geometry_normal = -hit_record.geometry_normal;
  if (dot(direction, hit_record.geometry_normal) * dot(direction, hit_record.normal) < 0.0) hit_record.normal = -hit_record.normal;
  
  if (!(abs(dot(hit_record.tangent, hit_record.tangent) - 1) < eps && abs(dot(hit_record.tangent, hit_record.normal)) < eps))
  {
    hit_record.tangent = gao1(hit_record.normal);
  }
  
  hit_record.tex_coord = mat3x2(v0.tex_coord, v1.tex_coord, v2.tex_coord) * ray_payload.barycentric;

  Material mat = materials[hit_record.hit_entity_id];

  if (mat.normal_texture_id >= 0) {
    hit_record.normal = normalize(mat3(hit_record.tangent,
        cross(hit_record.normal, hit_record.tangent), hit_record.normal) *
        (2.0f * vec3(texture(texture_samplers[mat.normal_texture_id], hit_record.tex_coord)) - vec3(1.0f)));
    hit_record.tangent = normalize(hit_record.tangent - dot(hit_record.tangent, hit_record.normal) * hit_record.normal);
  }

  return hit_record;
}