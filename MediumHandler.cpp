#include "core_settings.h"

//Prototype Volumetric Raytracer

HitInfo Raytracer::HandleRayMediumHit(Ray* ray, MediumIntersection mi, bool shadowRay)
{
	HitInfo hi = HitInfo();
	MediumBox mb = MediumBox();
	if (mi.mediumID >= 0)
		mb = *mediumBoxes[mi.mediumID];

	float mediumLength = length(mi.ls.end - mi.ls.start);
	//For a uniform chance r = [0,1], decide if this photon will interact with the medium, e.g. the attenuation distance
	float r = Rand(1.0f);

	//For light traveling in a medium;
	//float attenuation = exp(-extinction_coefficient * mediumLength);

	//pdf of this attenuation:
	//float pdf_attenunation = extinction_coefficient * exp(-extinction_coefficient * s);

	//integrating over the pdf from 0 to mediumlength to optain P(mediumLength), the cdf.
	//float P_attenuation = 1 - exp(-extinction_coefficient * s);

	//P_attenuation tells us the chance after travelling distance x of not having attenuated
	//Rewritting P_attenuation to isolate s, and inserting rnd gives us a distance of attenuation for a given probability;
	//float s = -log(rnd / extinction_coefficient);

	MediumParams mp;
	float extinction_coefficient;

	//Calculate s over the entire interval if homogenous medium, or perform ray-marching to find an attenuation distance
	//Get the corresponding medium parameters
	float s;
	if (mb.homogenous)
	{
		mp = mb.globalParams;
		extinction_coefficient = mp.absorption_coefficient + mp.scattering_coefficient;
		//Accumulate the extinction_coefficient of the medium 
		mediumExtinctionCoefficient[ray->pixelIndex].x += extinction_coefficient;
		s = -log(r) / extinction_coefficient;
	}
	else
	{
		float P_attenuation = 0;

		//Ray marching idea, accumulate P_attenuation - Pa, until we have crossed our goal attenuation r
		//Interpolate back to find the exact distance s_r where this attenuation is reached
		//Let T be (extinction_coefficient * s) in P_attenuation = 1 - exp(-extinction_coefficient * s) = 1 - exp(-T)
		//For the start of the march, P_attenuation_0 = 1 - exp(-T_0) = 0 => exp(-T_0) = 1; T_0 = 0;
		//For the goal of the march, P_attenuation_r = 1 - exp(-T_r) = r => exp(-T_r) = 1-r => T_r = -log(1-r);

		//We thus add to T_0 as we march the ray until we see T_0 pass T_r. Since T_0 = extinction_coefficient * s
		//Every march we add sample(extinction_coefficient,s) * mediumSampleStepSize to T_0 additionally
		//if we pass T_r with T_1 we find T_0 - T_r = extinction_coefficient * s - extinction_coefficient * s_r
		//=> (extinction_coefficient * s - extinction_coefficient * sGoal) / extinction_coefficient = s - s_r
		//So we find s_r by s -= s-s_r == s -= T_0 - T_r

		//To prevent patterns and allowing us to have a generally bigger step size so less sample, pick the 
		//extinction coefficient inside a step interval at a random location in that interval. 
		//So basically we do stratified monte carlo density sampling here.

		float T_0 = 0;
		float T_r = -log(r);
		float sampleLoc;
		float3 intervalLoc = ray->pos + mi.dist * ray->dir;
		for (s = 0; T_0 < T_r && s < mediumLength; s += mb.mediumSampleStepSize)
		{
			if (mb.stratifiedSampling)
				r = Rand(1.0f);
			else
				r = 1;
			mp = mb.SampleMedium(intervalLoc + r * mb.mediumSampleStepSize * ray->dir);
			intervalLoc += mb.mediumSampleStepSize * ray->dir;
			extinction_coefficient = mp.absorption_coefficient + mp.scattering_coefficient;
			//Accumulate the extinction_coefficient of the medium 
			mediumExtinctionCoefficient[ray->pixelIndex].x += extinction_coefficient;
			T_0 += mb.mediumSampleStepSize * extinction_coefficient;
		}
		//if(extinction_coefficient > 0)
		//	s -= (T_0 - T_r) / extinction_coefficient;

		//Know we know where we interact with the medium, sample this location for our final parameters
		mp = mb.SampleMedium(ray->pos + (mi.dist + s) * ray->dir);
	}

	//pdf of this attenuation:
	float pdf_attenunation = extinction_coefficient * exp(-extinction_coefficient * s);

	//Posteriori chance to reach the surface
	float P_surface = exp(-extinction_coefficient * s);

	//Now we have found for this photon where is will interact with the medium, assuming the medium is of infinite length
	//If the interaction point is actually outside of the medium, the photon does not interact and passes along.
	if (s > mediumLength)
	{
		//Apply emission
		ray->contribution += mp.emission_coefficient * mediumLength;
		ray->contribution;//The chance we 
		return hi;
	}

	//We are interacting with the medium
	hi.hitType = MEDIUM;

	//Shadow ray has interacted with the medium
	if (shadowRay)
		return hi;

	//The photon is guaranteed to be absorped, the ray dies here
	if (mp.scattering_coefficient == 0)
	{
		colorSampleData[ray->pixelIndex]++;
		return hi;
	}

	float scatterFrac = (mp.scattering_coefficient / extinction_coefficient);
	//Isotropic scattering; uniform random direction in 3D space
	if (mb.scatter_type == ISOTROPIC)
	{
		float rx = (Rand(1.0f) - 0.5) * 2;
		float ry = (Rand(1.0f) - 0.5) * 2;
		float rz = (Rand(1.0f) - 0.5) * 2;
		float mag = sqrt(rx * rx + ry * ry + rz * rz);
		float3 randomDir = float3{ rx / mag, ry / mag, rz / mag };

		//pdf_scattering = (mp.scattering_coefficient / extinction_coefficient) * phaseFunc
		//The phase function serves as a brdf but is also a pdf by definition, so division cancels out
		//For isotropic scattering, pdf = BRDF = phaseFunc = (1/(4*PI) so they cancel out, but well leave them in for clarity for now
		float3 brdf = ray->contribution * 1 / (4 * PI) * scatterFrac;
		float pdf_scattering = 1 / (4 * PI) * scatterFrac;

		ray->dist = mi.dist + s;
		ray->pos = ray->pos + ray->dist * ray->dir;
		ray->dir = randomDir;

		if (direct_sampling)
		{
			DirectSampleScatter(ray, scatterFrac, ISOTROPIC, 0);
		}

		//apply emission
		ray->contribution = brdf / pdf_scattering + mp.emission_coefficient * s;
	}
	else if (mb.scatter_type == ANISOTROPIC_HENYEY_GREENSTEIN)
	{
	https://www.astro.umd.edu/~jph/HG_note.pdf

	//Probability density function of Henyey_Greenstein for a direction theta
	//Given a modifier -1 <= g <= 1 ranging from backwards to forwards scattering  
	//p(theta) = (1/(4*PI)) * (1-g^2)/(1 + g^2 - 2g * cos(theta))^(-3/2)
	//Written as a function of u = cos(theta) =>
	//p(u) = (1/2) * (1-g^2)/(1 + g^2 - 2g * u)^(-3/2)

	//Accumulate the pdf so we have the cdf =>
	//P(u) = (1 - g^2)/2g *((1+g^2 - 2g * u)^(-1/2) - (1 + g)^-1))
	//u as a function of P =>
	//u = 1/2g*(1+g^2 - ((1-g^2)/1+gs) where s = 2P-1

	//We see that as P varies from 0-1, s varies from -1 - 1, and u ranges from -1 - 1.If
	//we then replace P by some r drawn uniformly at random on the interval[0, 1], the distribution
	//of the values of - will, for a large sample, approach the Henyey - Greenstein phase function.

		float r = Rand(1.0f);
		float g = mb.phaseFunctionParam;
		float p = 2 * r - 1;
		float u = (1 / (2 * g)) * (1 + g * g - pow((1 - g * g) / (1 + g * p), 2));	//Probably want to precompute this thing, once we start optimizing
		//Since u = cos(theta), theta = acos(u)
		float theta = acos(u);

		//Theta is now the angle the new ray is moving in wrt to the old direction, imagine the cone of directions around ray->dir for which the angle is theta
		//Assuming a polar coordinate system ,we thus take theta as the inclination, and we now need any phi for the azimuth to represent this cone of angles
		r = Rand(1.0f);
		float phi = r * 2 * PI;

		//Get the cartesian coordinates
		float3 nd = make_float3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));

		//Matrix multiplication to transform to world space
		float3 transformed = TangentToWorld(nd, ray->dir);

		//The phase function serves as a brdf but is also a pdf by definition, so division cancels out
		//so, no need to actually calculate the pdf of this function, so we wont in this case
		float3 brdf = ray->contribution * scatterFrac; //*phaseFunc
		float pdf_scattering = (mp.scattering_coefficient / extinction_coefficient); //*phaseFunc

		ray->dist = mi.dist + s;
		ray->pos = ray->pos + ray->dist * ray->dir;
		ray->dir = transformed;

		if (direct_sampling)
		{
			DirectSampleScatter(ray, scatterFrac, ANISOTROPIC_HENYEY_GREENSTEIN, g);
		}

		//apply emission
		ray->contribution = brdf / pdf_scattering + mp.emission_coefficient * s;
	}

	//Add scatter ray to the queue
	Ray nr;
	nr.dir = ray->dir;
	nr.pos = ray->pos + EPSILON * nr.dir;
	nr.contribution = ray->contribution;
	nr.pixelIndex = ray->pixelIndex;
	nr.recursionCount = ray->recursionCount + 1;
	nr.medium = ray->medium;
	nr.mat = ray->mat;
	nr.rayType = 2; //ScatterRayID
	AddRay(nr);

	return hi;
}

MediumParams MediumBox::SampleMedium(float3 pos)
{
	pos += mediumOffset;

	//Query the noise function, build the parameters
	float px = pos.x / zoom;
	float py = pos.y / zoom;
	float pz = pos.z / zoom;
	int tx = (int)(px) % textureSize;
	int ty = (int)(py) % textureSize;
	int tz = (int)(pz) % textureSize;

	if (tx < 0)
		tx *= -1;
	if (ty < 0)
		ty *= -1;
	if (tz < 0)
		tz *= -1;

	float density = noiseTexture[tx + textureSize * (ty + textureSize * tz)];
	//Average the density over the neighbouring noise pixels, may want to use guassian kernel
	//But considering the nature of noise this is probably kinda fine.
	if (trilinearInterpolation)
	{
		//Find the three noise cells that with the current cell form the cube in which our position lies, then interpolate the density 
		int h, v, d;

		//Get the distance along each direction
		float txd = px - (int)px;
		float tyd = py - (int)py;
		float tzd = pz - (int)pz;
		if (txd <= 0)
		{
			txd *= -1;
			h = (tx - 1) % textureSize;
		}
		else
			h = (tx + 1) % textureSize;
		if (tyd <= 0)
		{
			tyd *= -1;
			v = (ty - 1) % textureSize;
		}
		else
			v = (ty + 1) % textureSize;
		if (tzd <= 0)
		{
			tzd *= -1;
			d = (tz - 1) % textureSize;
		}
		else
			d = (tz + 1) % textureSize;

		if (h < 0)
			h *= -1;
		if (v < 0)
			v *= -1;
		if (d < 0)
			d *= -1;

		//Trilinear interpolation
		float v000_density = density;
		float v100_density = noiseTexture[h + textureSize * (ty + textureSize * tz)];
		float v010_density = noiseTexture[tx + textureSize * (v + textureSize * tz)];
		float v001_density = noiseTexture[tx + textureSize * (ty + textureSize * d)];
		float v101_density = noiseTexture[h + textureSize * (v + textureSize * tz)];
		float v011_density = noiseTexture[tx + textureSize * (v + textureSize * d)];
		float v110_density = noiseTexture[h + textureSize * (ty + textureSize * d)];
		float v111_density = noiseTexture[h + textureSize * (v + textureSize * d)];

		density = v000_density * (1 - txd) * (1 - tyd) * (1 - tzd) +
			v100_density * txd * (1 - tyd) * (1 - tzd) +
			v010_density * (1 - txd) * tyd * (1 - tzd) +
			v001_density * (1 - txd) * (1 - tyd) * tzd +
			v101_density * txd * (1 - tyd) * tzd +
			v011_density * (1 - txd) * tyd * tzd +
			v110_density * txd * tyd * (1 - tzd) +
			v111_density * txd * tyd * tzd;
		if (density < 0 || density > 1)
			cout << "bad density: " << density << " \n";
	}
	if (density < 0.25) //cutoff
		density = 0;

	MediumParams mp = MediumParams();
	mp.absorption_coefficient = globalParams.absorption_coefficient * density; //typical value for clouds
	mp.scattering_coefficient = globalParams.scattering_coefficient * density;
	mp.emission_coefficient = globalParams.emission_coefficient * density;
	return mp;
}

//NEE with a brdf for isotropicscatter -- reworking NEE to using pdf's prob gonna be nicer
void Raytracer::DirectSampleScatter(Ray* ray, float scatterFrac, ScatterType st = ISOTROPIC, float g = 0.5)
{
	int lights = triLights.size();
	int rndLight = RandomUInt() % lights;
	CoreLightTri* ctl = triLights[rndLight];
	float3 lp = RandomPointOnLight(ctl, ray);
	float3 iPos = ray->pos;
	float3 l = lp - iPos;
	float d = length(l);
	l = l / d;
	float cos_o = dot(-l, ctl->N);
	//Check if in correct direction
	if (cos_o <= 0) return;

	//ShadowRay
	Ray r = Ray();
	r.pos = iPos + EPSILON * l;
	r.dir = l;
	r.dist = d - 2 * EPSILON;


	HitInfo hit = RayIntersection(&r, true);
	if (hit.hitType != NONE) return;

	//Calc and set the color contributdddion of the next ray
	Material m = materials[ray->mat];
	float diffusion = 1 - m.specularity - m.glass - m.opacity;
	float solidAngle = cos_o * ctl->area / (d * d);
	float PDF = (1 / solidAngle) * scatterFrac; //t(x,z)
	float3 contribution;
	float BRDF;

	//The attenuation to reach the light source -> build into the original ray and shadow ray which has a probability of canceling the direct illumination
	//And the ray itself; Disadvantage: now we have done all this work for naught, if we can get an actual attenuation value that would be better
	//Advantage, No need to modify the shadow ray algo and fiddle with attenuation pdf, so leaving it be for now.
	if (st == ANISOTROPIC_HENYEY_GREENSTEIN)
	{
		//Calculate the phase function e.g. the brdf, this time the pdf is different however since we arent choosing a direction based on the phase function
		//Instead calculating the result of the phase function for a given direction theta with a certain PDF:
		//p(theta) = (1/(4*PI)) * (1-g^2)/(1 + g^2 - 2g * cos(theta))^(-3/2)

		//cos(theta) = dot(ray->dir, l) = u
		float u = dot(ray->dir, l);
		BRDF = (1 / (4 * PI)) * ((1 - g * g) / pow((1 + g * g - 2 * g * u), (-3 / 2))) * scatterFrac;
		contribution = ray->contribution * BRDF * ctl->radiance * diffusion * lights / PDF;
	}
	else
	{
		//Again, scatterFrac actually cancels out here, but leaving it in for clarity for now
		BRDF = (1 / (4 * PI)) * scatterFrac;
		contribution = ray->contribution * BRDF * ctl->radiance * diffusion * lights / PDF;
	}
	ApplyContribution(ray, contribution, MEDIUM_ILLUMINATION);
}

void MediumBox::BuildNoiseTexture()
{
	noiseTexture = new float[textureSize * textureSize * textureSize];

	cout << "Building Noise Map... \n";
	cout << "Building Worley Noise... \n";
	//First layer of Cellular noise
	FastNoiseLite noise;
	noise.SetNoiseType(FastNoiseLite::NoiseType_Cellular);
	float n;
	float zf = 1; //zoomFactor
	for (int x = 0; x < textureSize; x++)
		for (int y = 0; y < textureSize; y++)
			for (int z = 0; z < textureSize; z++)
			{
				n = (-noise.GetNoise(x * zf, y * zf, z * zf) + 1) / 2;
				noiseTexture[x + textureSize * (y + textureSize * z)] = n;
			}

	float blendFactor = 0.5; //Averaging factor, contribution of new noise
	zf = 1;

	//Second layer of Perlin noise
	cout << "Building Perlin noise\n";
	noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
	noise.SetFrequency(0.01);
	noise.SetFractalType(FastNoiseLite::FractalType_FBm);
	noise.SetFractalGain(0.7);
	noise.SetFractalOctaves(10);
	for (int x = 0; x < textureSize; x++)
		for (int y = 0; y < textureSize; y++)
			for (int z = 0; z < textureSize; z++)
			{
				n = (-noise.GetNoise(x * zf, y * zf, z * zf) + 1) / 2;
				noiseTexture[x + textureSize * (y + textureSize * z)] = noiseTexture[x + textureSize * (y + textureSize * z)] * (1 - blendFactor) * blendFactor * n;
			}

	//Contrast adjustment
	float maxVal = 0;
	for (int x = 0; x < textureSize; x++)
		for (int y = 0; y < textureSize; y++)
			for (int z = 0; z < textureSize; z++)
			{
				float newVal = pow(noiseTexture[x + textureSize * (y + textureSize * z)], 2);
				if (newVal > maxVal)
					maxVal = newVal;
				noiseTexture[x + textureSize * (y + textureSize * z)] = newVal;
			}

	//Normalize and threshold low points
	for (int x = 0; x < textureSize; x++)
		for (int y = 0; y < textureSize; y++)
			for (int z = 0; z < textureSize; z++)
			{
				float newVal = noiseTexture[x + textureSize * (y + textureSize * z)] / maxVal;
				noiseTexture[x + textureSize * (y + textureSize * z)] = newVal;
				//Threshold the noise to give the texture sharp edges and create bigger gaps
				if (newVal < 0.50)
				{
					//newVal = pow(newVal, 2);
					noiseTexture[x + textureSize * (y + textureSize * z)] = 0;
				}
			}
	cout << "Done. \n";
}

//Returns the segment of the ray intersecting the aabb medium box
MediumIntersection MediumBox::Intersect(Ray* r)
{
	MediumIntersection mi = MediumIntersection();
	float3 o = start;
	float3 e = end;

	float tmin = (o.x - r->pos.x) / r->dir.x;
	float tmax = (e.x - r->pos.x) / r->dir.x;

	if (tmin > tmax)
	{
		swap(tmin, tmax);
	}

	float tymin = (o.y - r->pos.y) / r->dir.y;
	float tymax = (e.y - r->pos.y) / r->dir.y;

	if (tymin > tymax)
	{
		swap(tymin, tymax);
	}

	if ((tmin > tymax) || (tymin > tmax))
	{
		return mi;
	}

	if (tymin > tmin)
	{
		tmin = tymin;
	}

	if (tymax < tmax)
	{
		tmax = tymax;
	}

	float tzmin = (o.z - r->pos.z) / r->dir.z;
	float tzmax = (e.z - r->pos.z) / r->dir.z;

	if (tzmin > tzmax)
	{
		swap(tzmin, tzmax);
	}

	if ((tmin > tzmax) || (tzmin > tmax))
	{
		return mi;
	}

	if (tzmin > tmin)
	{
		tmin = tzmin;
	}

	if (tzmax < tmax)
	{
		tmax = tzmax;
	}

	if (tmin <= 0 && tmax <= 0)
		return mi;

	LineSegment ls = LineSegment();
	ls.start = r->pos + r->dir * tmin;
	ls.end = r->pos + r->dir * tmax;

	mi.ls = ls;
	mi.dist = tmin;
	return mi;
}

//Does not support other objects intersecting with this box
MediumIntersection Raytracer::GetNearestMediumBoxIntersection(Ray* ray)
{
	MediumIntersection mi = MediumIntersection();
	for (int i = 0; i < mediumBoxes.size(); i++)
	{
		MediumIntersection mediumIntersection = mediumBoxes[i]->Intersect(ray);
		if (mediumIntersection.dist < ray->dist)
		{
			mi.ls = mediumIntersection.ls;
			mi.dist = mediumIntersection.dist;
			mi.mediumID = i;
		}
	}
	return mi;
}