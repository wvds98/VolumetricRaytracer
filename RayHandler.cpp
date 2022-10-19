#include "core_settings.h"

inline void async_trace(int id, Raytracer* raytracer, Ray r)
{
	raytracer->Trace(&r);
}

Ray Raytracer::InitRaySphere(const ViewPyramid& view, float x, float y, int index, float3 forward, float3 left, float3 down)
{
	float4 middle = normalize(make_float4(forward, 0));
	mat4 rotation_y = mat4::Rotate(down, fish_eye_x_angle * x);
	mat4 rotation_x = mat4::Rotate(left, fish_eye_y_angle * y);
	float4 dir = (middle * rotation_y) * rotation_x;

	Ray r;

	r.pos = view.pos;
	r.dir = make_float3(dir);
	r.dist = INFINITY;

	//Set the pixel this ray belongs to
	r.pixelIndex = index;
	r.contribution = make_float3(1);
	return r;
}

Ray Raytracer::InitRay(const ViewPyramid& view, float x, float y, int index, float3 forward, float3 left, float3 down)
{
	////Calculate the ray direction
	float3 p = forward + x * left + y * down;
	float3 d = normalize(p);

	//Set ray origin and position
	Ray r;

	r.pos = view.pos;
	r.dir = d;
	r.dist = INFINITY;

	//Set the pixel this ray belongs to
	r.pixelIndex = index;
	r.contribution = make_float3(1);
	return r;
}

float3 Raytracer::GetTriangleIntersection(float4 v1, float4 v2, float4 v3, Ray* ray)
{
	const float epsilon = 0.0000001;
	float3 e1 = make_float3(v2 - v1);
	float3 e2 = make_float3(v3 - v1);
	//Find intersection point
	float3 h = cross(ray->dir, e2);
	float a = dot(e1, h);

	//Check parallel case
	if (a > -epsilon && a < epsilon)
	{
		return make_float3(-1);
	}

	float f = 1.0 / a;
	float3 s = ray->pos - make_float3(v1);
	float u = f * dot(s, h);

	if (u < 0.0 || u > 1.0)
	{
		return make_float3(-1);
	}

	float3 q = cross(s, e1);
	float v = f * dot(ray->dir, q);

	if (v < 0.0 || u + v > 1.0)
	{
		return make_float3(-1);
	}

	float t = f * dot(e2, q);

	return make_float3(t, u, v);
}

HitInfo Raytracer::RayIntersection(Ray* ray, bool shadowRay = false)
{
	float maxDist = ray->dist;
	float3 baseContribution = ray->contribution;
	bool hit = false;

	if (intersect_primitives)
		hit = GetNearestPrimitiveIntersection(primitives, ray, baseContribution, shadowRay);

	//Only intersect the TriLights if pathtracing and not using 
	//direct sampling or if the ray is a primary ray or a direct ray.
	CoreLightTri* ctl = nullptr;
	if (path_trace && !hit && (ray->rayType == 1 || ray->recursionCount == 0 || !direct_sampling || MIS) && !shadowRay)
		ctl = GetNearestTriLightIntersection(ray);

	if (intersect_meshes)
	{
		if(useBVH)
			hit = GetNearestBVHIntersection(ray, baseContribution, shadowRay);
		else
			hit = GetNearestMeshIntersection(ray, baseContribution, shadowRay);
	}

	//Check if we interact with a medium before hitting the found closest object
	//Handle the medium boxes -- Boxes are not allowed to overlap for now for simplicity, a ray can only intersect 1 medium
	if (intersect_mediums && ray->medium == 0)//Check the ray is not travelling through a solid medium
	{
		MediumIntersection mi;
		if(!globalMedium)
			mi = GetNearestMediumBoxIntersection(ray);
		else
		{
			mi.ls.start = ray->pos;
			mi.ls.end = ray->pos + ray->dir * ray->dist;
			mi.dist = 0;
			mi.mediumID = -2;
		}
		if (mi.mediumID != -1)
		{			HitInfo hi = HandleRayMediumHit(ray,mi, shadowRay);
			if (hi.hitType == MEDIUM)
			{
				return hi;
			}
		}
	}

	//Check if we hit a light source first
	//If so apply it's color and stop handling this ray
	if (ctl != nullptr && !hit && path_trace)
	{
		ray->contribution = baseContribution;
		HitInfo hi = HitInfo();
		hi.clt = ctl;
		hi.hitType = LIGHT;
		return hi;
	}

	if (ray->dist < maxDist)
	{
		HitInfo hi = HitInfo();
		hi.hitType = OBJECT;
		return hi;
	}
	else
	{
		return HitInfo();
	}
}

void Raytracer::GetSkyDome(Ray* ray)
{
	float3 tDir = -worldToLight.TransformVector(ray->dir);
	float3 skyPixel = SampleSkyDome(tDir); // optionally select the scaled skydome for efficiency here
	float3 contribution = skyPixel * ray->contribution * skyDomeBrightNess;
	//Store primary hits in albedo map
	if (ray->recursionCount == 0 && splitAlbedo)
	{
		albedoData[ray->pixelIndex] += contribution;
		contribution = contribution;
	}
	ApplyContribution(ray, contribution);
	colorSampleData[ray->pixelIndex]++;
}

void Raytracer::ApplyTriLight(Ray* ray, CoreLightTri* ctl)
{
	auto radiance = ctl->radiance;
	if (!path_trace)
		radiance *= ctl->energy / (ray->dist * ray->dist);

	//Calculate and apply the color contribution of this ray to the pixel.
	if (MIS && ray->rayType == 0 && ray->recursionCount != 0) //MIS cant be used after a specular vertex or scatter rays, rayType must be 0, and cant be a primary ray;
	{
		//When this ray was created the brdf was added to it at full strength
		float brdfPDF = ray->prevPDF;

		//Calculate lightPDF
		float3 lp = ray->pos + ray->dist * ray->dir;
		float3 iPos = ray->pos;
		float3 l = lp - iPos;
		float d = length(l);
		l = l / d;
		float cos_o = dot(-l, ctl->N);
		float cos_i = dot(l, ray->norm);

		float solidAngle = cos_o * ctl->area / (d * d);
		float lightPDF = 1 / solidAngle;

		//ray->contribution = randomWalkContribution / brdfPDF
		//must become; randomWalkContribution / (brdfPDF + lightPDF)
		//randomWalkContribtuion / (brdfPDF + lightPDF)
		float3 randomWalkContribution = ray->contribution * brdfPDF;
		float PDF = brdfPDF + lightPDF;
		float3 newContribution = randomWalkContribution / PDF;
		ray->contribution = newContribution;
	}
	if (ray->recursionCount == 0 && splitAlbedo)
	{
		RecordPrimaryHit(ray, ray->contribution * radiance);
	}
	ApplyContribution(ray, ray->contribution * radiance);
	colorSampleData[ray->pixelIndex]++;
}

void Raytracer::Trace(Ray* ray)
{
	HitInfo hit = RayIntersection(ray);

	if (hit.hitType == NONE)
	{
		GetSkyDome(ray);
	}
	else if(hit.hitType == OBJECT) 
	{
		HandleRayHit(ray);
	}
	else if (hit.hitType == LIGHT) 
	{
		ApplyTriLight(ray, hit.clt);
	}
}

bool Raytracer::GetNearestPrimitiveIntersection(vector<Primitive*> primitives, Ray* ray, float3 baseContribution, bool shadowRay = false)
{
	Primitive* closest_primitive = nullptr;
	for (auto primitive : primitives)
	{
		float t = primitive->GetIntersectionDistance(ray);
		if (t > 0)
		{
			float closest = ray->dist;
			if (t < closest)
			{
				ray->dist = t;
				closest_primitive = primitive;

				if (shadowRay) //Any valid intersection, stop looking for intersections
				{
					//Check if the intersected object is not glass
					//For now, shadows pass straight through glass and objects with opacity
					Material m = materials[closest_primitive->mat];
					if (m.glass == 0 && m.opacity == 0)
					{
						break;
					}
					else
					{
						//Forget this intersection
						ray->dist = closest;
						closest_primitive = nullptr;
					}
				}
			}
		}
	}
	if (closest_primitive != nullptr)
	{
		auto pos = ray->pos + ray->dir * ray->dist;
		ray->norm = closest_primitive->GetNormal(pos);
		ray->mat = closest_primitive->mat;

		//Temp code hack for applying a pattern to plane until we get textures
		Material mat = materials[ray->mat];
		float3 color = materials[ray->mat].color;
		float3 iPos = ray->pos + ray->dist * ray->dir;
		if (mat.specularity > 0.2 && mat.specularity < 0.21)
		{
			bool black = (((int)(iPos.x / 2) % 2) + ((int)(iPos.z / 2) % 2)) % 2 == 0;
			if ((iPos.x < 0 || iPos.z < 0) && (iPos.x > 0 || iPos.z > 0))
			{
				black = !black;
			}
			if (black)
			{
				color = { 0.3,0.3,0.3 };
			}
			else
			{
				color = { 1,1,1 };
			}
		}

		ray->contribution = baseContribution * color;

		return true;
	}
	return false;
}

CoreLightTri* Raytracer::GetNearestTriLightIntersection(Ray* ray)
{
	CoreLightTri* closest_light = nullptr;
	for (int i = 0; i < triLights.size(); i++)
	{
		CoreLightTri* tl = triLights[i];
		float3 ti = GetTriangleIntersection(make_float4(tl->vertex0), make_float4(tl->vertex1), make_float4(tl->vertex2), ray);
		if (ti.x > 0)
		{
			float closest = ray->dist;
			if (ti.x < closest)
			{
				ray->dist = ti.x;
				closest_light = triLights[i];
			}
		}
	}
	return closest_light;
}

bool Raytracer::GetNearestBVHIntersection(Ray* ray, float3 baseContribution, bool shadowRay = false)
{
	int ct = -1;
	CoreTri* closest_triangle = nullptr;
	float closest_u = 0;
	float closest_v = 0;
	for (int i = 0; i < meshes.size(); i++)
	{
		Mesh m = meshes[i];
		BVHInstance bvhi = bvhInstances[i];
		bvhi.bvh->TraverseBVH(ray, bvhi.bvh->root, &bvhi.transform, &m, shadowRay, &ct, &closest_u, &closest_v, &materials, baseMaterialCount);
		if (ct != -1)
			closest_triangle = &m.triangles[ct];
	}

	//Process triangle data
	if (closest_triangle != nullptr)
	{
		//We don't want to cull triangle faces, so check if this normal lies in 
		//the direction of the ray, if it does flip the normal, unless we are inside the mesh
		ray->norm = { closest_triangle->Nx, closest_triangle->Ny, closest_triangle->Nz };
		if (dot(ray->norm, ray->dir) > 0 && ray->medium == 0)
		{
			ray->norm *= -1;
		}

		//Fetch material
		ray->mat = closest_triangle->material + baseMaterialCount;

		ReadMaterial(ray, closest_triangle, baseContribution, closest_u, closest_v);

		return true;
	}
	return false;
}

bool Raytracer::GetNearestMeshIntersection(Ray* ray, float3 baseContribution, bool shadowRay = false)
{
	CoreTri* closest_triangle = nullptr;
	float closest_u = 0;
	float closest_v = 0;
	for (int i = 0; i < meshes.size(); i++) for (int j = 0; j < meshes[i].vcount; j += 3) //meshes[i].vcount
	{
		Mesh m = meshes[i];

		float4 v1 = m.vertices[j];
		float4 v2 = m.vertices[j + 1];
		float4 v3 = m.vertices[j + 2];
		float3 ti = GetTriangleIntersection(v1, v2, v3, ray);
		if (ti.x > 0)
		{
			float closest = ray->dist;
			if (ti.x < closest)
			{
				ray->dist = ti.x;
				closest_u = ti.y;
				closest_v = ti.z;
				closest_triangle = &meshes[i].triangles[j / 3];

				Material m = materials[closest_triangle->material + baseMaterialCount];

				//Skip emmissive triangles (They shouldnt be inside meshes) -- temporary solution
				if (m.cmat->color.value.x > 1 || m.cmat->color.value.y > 1 || m.cmat->color.value.z > 1)
				{
					//Forget this intersection
					ray->dist = closest;
					closest_triangle = nullptr;
					continue;
				}

				if (shadowRay) //Any valid intersection, stop looking for intersections
				{
					//Check if the intersected object is not glass
					//For now, shadows pass straight through glass, opacity and emmision. 				
					if (m.glass == 0 && m.opacity == 0)
					{
						break;
					}
					else
					{
						//Forget this intersection
						ray->dist = closest;
						closest_triangle = nullptr;
					}
				}
			}
		}
	}

	//Process triangle data
	if (closest_triangle != nullptr)
	{
		//We don't want to cull triangle faces, so check if this normal lies in 
		//the direction of the ray, if it does flip the normal, unless we are inside the mesh
		ray->norm = { closest_triangle->Nx, closest_triangle->Ny, closest_triangle->Nz };
		if (dot(ray->norm, ray->dir) > 0 && ray->medium == 0)
		{
			ray->norm *= -1;
		}

		//Fetch material
		ray->mat = closest_triangle->material + baseMaterialCount;

		ReadMaterial(ray, closest_triangle, baseContribution, closest_u, closest_v);

		return true;
	}
	return false;
}

void Raytracer::HandleRayHit(Ray* ray)
{
	float reflectionCutoff = 0.01f;
	float glassinessCutoff = 0.01f;
	Material mat = materials[ray->mat];
	Material medium = materials[ray->medium];
	float specularity = mat.specularity;
	bool glass = mat.glass;
	float opacity = mat.opacity;

	//If this ray doesn't travel through air
	//Apply absorption if applicable
	if (medium.refractionIndice != 1.0)
	{
		float3 a = medium.absorption;
		float3 I = ray->contribution;
		float d = ray->dist;
		if (a.x != 0 || a.y != 0 || a.z != 0)
		{
			I.x *= exp(-a.x * d);
			I.y *= exp(-a.y * d);
			I.z *= exp(-a.z * d);
			ray->contribution = I;
		}
	}

	//If path tracing, take a single path based on probability dist
	//Based on assumption specularity + opacity + diffuse + glass = 1;
	//Special case, opacity applied to glass is legal
	if (path_trace)
	{
		float rnd = Rand(1.0f);
		if (rnd < specularity)
		{
			HandleReflection(ray, 1);
		}
		else if (glass)
		{
			HandleGlassIntersection(ray);
		}
		else if (rnd > specularity && rnd <= specularity + opacity)
		{
			HandleOpacity(ray, 1);
		}
		else if (rnd > specularity + opacity)
		{
			HandleDiffuse(ray);
		}
	}
	else
	{
		if (opacity > 0)
		{
			float opacity = materials[ray->mat].opacity;
			HandleOpacity(ray, opacity);
			ray->contribution *= 1 - opacity;
		}
		if (glass)
			HandleGlassIntersection(ray);
		else
		{
			if (specularity > reflectionCutoff && !glass)
				HandleReflection(ray, specularity);
			if (specularity < 1 - reflectionCutoff && !glass)
				HandleDiffuse(ray);
		}
	}

}

void Raytracer::HandleReflection(Ray* ray, float f)
{
	//ReflectVector
	Ray nr = ReflectRay(ray);
	nr.contribution = ray->contribution * f;
	nr.pixelIndex = ray->pixelIndex;
	nr.recursionCount = ray->recursionCount + 1;
	nr.rayType = 1;

	nr.medium = ray->medium;
	nr.mat = ray->mat;

	//nr.medium.refractionIndice = ray->medium.refractionIndice;
	//nr.medium.absorption = ray->medium.absorption;
	//nr.mat.cmat = ray->mat.cmat;

	AddRay(nr);
}

void Raytracer::HandleDiffuse(Ray* ray)
{
	if (!path_trace)
	{
		ApplyDirectIllumination(ray);
	}
	else
	{
		//Record albedo for primary rays
		if (ray->recursionCount == 0 && splitAlbedo)
		{
			RecordPrimaryHit(ray, ray->contribution);
			ray->contribution = { 1,1,1 };
		}

		if (direct_sampling)
		{
			DirectSample(ray);
		}

		if (russian_roulette)
		{
			//Russian roulette on throughput value, the lower the contribution, the less interesting
			//following this path becomes.
			float max1 = max(ray->contribution.x, ray->contribution.y);
			float max2 = max(max1, ray->contribution.z);
			float p_survive = clamp(max2, 0.1f, 0.9f);
			if (Rand(1.0f) > p_survive)
			{
				colorSampleData[ray->pixelIndex]++;
				return;
			}
			else
				ray->contribution /= p_survive;
		}

		if (!use_importance_sampled_randomwalk)
			RandomWalk(ray);
		else
			RandomWalkCosineWeightedDiffuseReflection(ray);
	}
}

void Raytracer::DirectSample(Ray* ray)
{
	int lights = triLights.size();
	int rndLight = RandomUInt() % lights;
	CoreLightTri* ctl = triLights[rndLight];
	float3 lp = RandomPointOnLight(ctl,ray);
	float3 iPos = ray->pos + ray->dir * ray->dist;
	float3 l = lp - iPos;
	float d = length(l);
	l = l / d;
	float cos_o = dot(-l, ctl->N);
	float cos_i = dot(l, ray->norm);
	//Check if in correct direction
	if (cos_o <= 0 || cos_i <= 0) return;

	//ShadowRay
	Ray r = Ray();
	r.pos = ray->pos + ray->dir * ray->dist + EPSILON * l;
	r.dir = l;
	r.dist = d - 2 * EPSILON;

	HitInfo hit = RayIntersection(&r, true);
	if (hit.hitType != NONE) return;

	//Calc and set the color contributdddion of the next ray
	Material m = materials[ray->mat];
	float diffusion = 1 - m.specularity - m.glass - m.opacity;

	float3 BRDF = ray->contribution * INVPI;			//why not 1/(2*PI)??
	float solidAngle = cos_o * ctl->area / (d * d);
	float PDF = 1 / solidAngle;
	if (MIS)
	{
		if (use_importance_sampled_randomwalk)
			PDF += cos_i / PI;
		else
			PDF += 1 / (2 * PI);
	}
	float3 contribution = BRDF * lights * ctl->radiance * cos_i * diffusion / PDF;
	ApplyContribution(ray, contribution);
}

void Raytracer::ApplyContribution(Ray* ray, float3 contribution, IlluminationType it)
{
	colorData[ray->pixelIndex] += contribution;

	if (it == DIFFUSE_ILLUMINATION)
	{
		illuminationData[ray->pixelIndex] += contribution;
	}
	else if (it == MEDIUM_ILLUMINATION)
	{
		mediumData[ray->pixelIndex] += contribution;
	}
}

void Raytracer::RecordPrimaryHit(Ray* ray, float3 albedo)
{
	albedoData[ray->pixelIndex] += albedo;
	normalData[ray->pixelIndex] += ray->norm;
	worldPosData[ray->pixelIndex] += ray->pos + ray->dist * ray->dir;
}

//Returns a uniform random point within the triangleb
//According to https://www.cs.princeton.edu/~funk/tog02.pdf 4.2
float3 Raytracer::RandomPointOnLight(CoreLightTri* ctl, Ray* ray)
{
	float r1 = Rand(1.0);
	float r2 = Rand(1.0);
	if (blueNoise != nullptr)
	{
		const uint x = (ray->pixelIndex % pixelData->width) & 127, y = ((ray->pixelIndex) / pixelData->width) & 127;
		float4 r4 = BlueNoiseSampler4(x, y, colorSampleData[ray->pixelIndex], 4 * ray->recursionCount - 4, blueNoise);
		r1 = r4.x;
		r2 = r4.y;
	}
	float sqrtr1 = sqrt(r1);
	float3 p = (1 - sqrtr1) * ctl->vertex0 + (sqrtr1 * (1 - r2)) * ctl->vertex1 + (r2 * sqrtr1) * ctl->vertex2;
	return p;
}

void Raytracer::RandomWalk(Ray* ray)
{
	//Get spherical coordinates of normal vector
	float theta = atan2(ray->norm.y, ray->norm.x);
	float phi = acos(ray->norm.z);

	//Set new spherical coordinates randomly over signed 90 degrees
	float ntheta = theta + ((0.5 - Rand(1.0f)) * PI);
	float nphi = phi + ((0.5 - Rand(1.0f)) * PI);

	//Calculate cartesian coordinates from new spherical coordinates
	float3 newDir = make_float3(sin(nphi) * cos(ntheta), sin(nphi) * sin(ntheta), cos(nphi));

	Ray nr;
	nr.dir = newDir;
	nr.pos = ray->pos + ray->dist * ray->dir + EPSILON * ray->norm;

	//Apply BRDF to find the new contribution of the ray
	//Diffusie BRDF component
	Material m = materials[ray->mat];
	float diffusion = 1 - m.specularity - m.glass - m.opacity;

	nr.pixelIndex = ray->pixelIndex;
	nr.recursionCount = ray->recursionCount + 1;
	nr.medium = ray->medium;
	nr.mat = ray->mat;

	//Calc and set the color contribution of the next ray
	float3 BRDF = ray->contribution * INVPI;
	float PDF = 1 / (2 * PI);
	nr.prevPDF = PDF;
	nr.contribution = BRDF * (dot(ray->norm, newDir)/PDF) * diffusion;
	AddRay(nr);
}

float3 Raytracer::TangentToWorld(float3 td, float3 n)
{
	//Generate tangent and bitangent
	float3 w = { 0,1,0 };
	if (n.y > 0.99f)
		w = { 1,0,0 };
	float3 t = cross(n, w);
	float3 b;
	if (length(t) == 0)
		t = cross(n, { 0,0,1 });
	t = normalize(t);
	b = normalize(cross(n, t));

	//Matrix multiplication to transform from tangent to world space
	float3 transformed = make_float3(t.x * td.x + b.x * td.y + n.x * td.z,
		t.y * td.x + b.y * td.y + n.y * td.z,
		t.z * td.x + b.z * td.y + n.z * td.z);

	return transformed;
}

//Importance sampled diffuse reflection for the random walk
//Code copied from slides source global illumination compendium
void Raytracer::RandomWalkCosineWeightedDiffuseReflection(Ray* ray)
{
	//Get cosine weighted direction
	float r0 = Rand(1.0);
	float r1 = Rand(1.0);
	float r = sqrt(r0);
	float theta = 2 * PI * r1;
	float x = r * cosf(theta);
	float y = r * sinf(theta);
	float3 cd = make_float3(x, y, sqrt(1 - r0)); 

	float3 transformed = TangentToWorld(cd, ray->norm);

	Ray nr;
	nr.dir = transformed;
	nr.pos = ray->pos + ray->dist * ray->dir + EPSILON * ray->norm;

	//Apply BRDF to find the new contribution of the ray
	//Diffusie BRDF component
	Material m = materials[ray->mat];
	float diffusion = 1 - m.specularity - m.glass - m.opacity;

	nr.pixelIndex = ray->pixelIndex;
	nr.recursionCount = ray->recursionCount + 1;
	nr.medium = ray->medium;
	nr.mat = ray->mat;

	//Calc and set the color contribution of the next ray
	float3 BRDF = ray->contribution * INVPI;
	float PDF = dot(ray->norm, transformed) * INVPI;
	nr.contribution = BRDF * dot(ray->norm, transformed) * diffusion / PDF;

	AddRay(nr);
}

//pass the light straight through
void Raytracer::HandleOpacity(Ray* ray, float f)
{
	//Copy ray
	Ray nr = Ray();
	nr.dir = ray->dir;
	float3 iPos = ray->pos + ray->dist * ray->dir;
	nr.pos = iPos + EPSILON * nr.dir;
	nr.contribution = ray->contribution * f;
	nr.pixelIndex = ray->pixelIndex;
	nr.recursionCount = ray->recursionCount + 1;

	nr.medium = ray->medium;
	nr.mat = ray->mat;

	//nr.medium.refractionIndice = ray->medium.refractionIndice;
	//nr.medium.absorption = ray->medium.absorption;
	//nr.mat.cmat = ray->mat.cmat;
	AddRay(nr);
}

void Raytracer::HandleGlassIntersection(Ray* ray)
{
	float n1 = materials[ray->medium].refractionIndice;

	float matr = materials[ray->mat].refractionIndice;
	//If we are inside an object, i.e, mediumRefrac != 1 (air), we assume the exterior is air
	//; n2 = 1.0. We also flip the normal of the intersection point
	float n2 = 1.0;
	if (n1 == 1.0)
	{
		n2 = matr;
	}
	else
	{
		ray->norm *= -1;
	}

	//Not sure why or when this can happen, but it shouldnt I don't think
	if (n1 == n2)
	{
		return;
	}

	float cos1 = dot(-ray->dir, ray->norm);
	float n = (n1 / n2);
	float k = 1 - n * n * (1 - cos1 * cos1);

	float oi = acos(cos1);
	float ia = matr * sin(oi);
	float ot = asin(ia / matr);

	//Fresnell Euqation possibly not correctly implemented
	float cos0t = sqrt(1 - (n * sin(oi)) * (n * sin(oi)));
	float eq1 = (n1 * cos1 - n2 * cos0t) / (n1 * cos1 + n2 * cos0t);
	float eq2 = (n1 * cos0t - n2 * cos1) / (n1 * cos0t + n2 * cos1);
	float fr = 0.5 * (eq1 * eq1 + eq2 * eq2);

	if (k < 0)
	{
		fr = 1;
	}

	//This can apparently happen, but it shouldnt, gotta look into that, probably
	//same reason n1 can be n2
	if (fr < 0 || fr > 1)
	{
		fr = clamp(fr, 0.0f, 1.0f);
	}

	//If path tracing, take a single path
	if (path_trace)
	{
		float rnd = Rand(1.0f);
		if (rnd < fr)
		{
			HandleReflection(ray, 1);
		}
		else
		{
			HandleTransmission(ray, n, n1, cos1, k, 0);
		}
	}
	else
	{
		HandleReflection(ray, fr);
		HandleTransmission(ray, n, n1, cos1, k, fr);
	}
	//Handle reflection component
}

void Raytracer::HandleTransmission(Ray* ray, float n, float n1, float cos1, float k, float f)
{

	//Handle transmission component
	Ray nr = Ray();
	float3 T = n * ray->dir + ray->norm * (n * cos1 - sqrt(k));
	nr.dir = T;
	float3 iPos = ray->pos + ray->dist * ray->dir;
	nr.pos = iPos + EPSILON * nr.dir;
	nr.contribution = ray->contribution * (1 - f);
	nr.pixelIndex = ray->pixelIndex;
	nr.recursionCount = ray->recursionCount + 1;
	nr.rayType = 1;

	nr.mat = ray->mat;
	//If ray exits a object, set medium material to default, which has air refraction: 0
	//and flip normal, on enter, set medium to material of the intersected object
	if (n1 != 1.0)
	{
		nr.medium = 0;
	}
	else
	{
		nr.medium = ray->mat;
	}

	AddRay(nr);
}

//Add the ray to the ray stack
//Can also call Trace on this ray instead
//In the case of a recursive implementation
void Raytracer::AddRay(Ray nr)
{
	if (nr.recursionCount < MAXRECURSION)
	{
		//rayStack.add(nr);
		thread_pool.push(async_trace, this, nr);
	}
	else if (path_trace)
	{
		//colorSampleData[nr.pixelIndex]++;
	}
}

//Calculate a new reflected ray.
Ray Raytracer::ReflectRay(Ray* ray)
{
	Ray r;
	r.dir = ray->dir - 2 * (dot(ray->dir, ray->norm)) * ray->norm;
	float3 iPos = ray->pos + ray->dist * ray->dir;
	r.pos = iPos + EPSILON * r.dir;
	return r;
}

void Raytracer::ApplyDirectIllumination(Ray* ray)
{
	float3 intersection_point = ray->pos + ray->dist * ray->dir;
	for (CorePointLight* light : pointLights)
	{
		auto direction = light->position - intersection_point;
		auto distance = length(direction);
		direction = normalize(direction);
		Ray shadowRay;
		shadowRay.dir = direction;
		shadowRay.dist = distance;
		shadowRay.pos = intersection_point + EPSILON * ray->norm;
		HitInfo hit = RayIntersection(&shadowRay, true);

		if (hit.hitType == NONE)
		{
			ApplyLightContribution(ray, light->radiance, direction, light->energy, distance);
		}
	}
	for (CoreDirectionalLight* light : directionalLights)
	{
		auto direction = -light->direction;
		Ray shadowRay;
		shadowRay.dir = direction;
		shadowRay.pos = intersection_point + EPSILON * ray->norm;
		HitInfo hit = RayIntersection(&shadowRay, true);

		if (hit.hitType == NONE)
		{
			ApplyLightContribution(ray, light->radiance, direction, light->energy, 1);
		}
	}
	for (CoreSpotLight* light : spotLights)
	{
		auto direction = light->position - intersection_point;
		auto distance = length(direction);
		direction = normalize(direction);
		Ray shadowRay;
		shadowRay.dir = direction;
		shadowRay.dist = distance;
		shadowRay.pos = intersection_point + EPSILON * ray->norm;
		HitInfo hit = RayIntersection(&shadowRay, true);

		if (hit.hitType == NONE)
		{
			//Calc if within the spotlights range
			float d = dot(-direction, light->direction);
			if (d <= light->cosOuter)
				break;
			float fallOffInterval = light->cosInner - light->cosOuter;
			float intervalFrac = (clamp(d, light->cosOuter, light->cosInner) - light->cosOuter);
			float fallOff = intervalFrac / fallOffInterval;

			float nDotL = dot(ray->norm, direction);
			float3 lr = light->radiance;
			auto energy = lr.x + lr.y + lr.z;
			ApplyLightContribution(ray, lr, direction, energy * fallOff, distance);
		}
	}
	colorSampleData[ray->pixelIndex]++;
}

void Raytracer::ApplyLightContribution(Ray* ray, float3 r, float3 direction, float energy, float distance)
{
	float nDotL = dot(ray->norm, direction);
	auto radiance = r * energy / (distance * distance) * nDotL;
	Material m = materials[ray->mat];

	//Calculate and apply the color contribution of this ray to the pixel.
	float3 contribution = ray->contribution * radiance *
		(1 - m.specularity) * (1 - m.glass) * (1 - m.opacity);
	ApplyContribution(ray, contribution);
}