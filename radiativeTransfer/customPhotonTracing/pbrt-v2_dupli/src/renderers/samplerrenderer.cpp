
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


// renderers/samplerrenderer.cpp*
#include "stdafx.h"
#include "renderers/samplerrenderer.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"

extern bool PhotonImage;

static uint32_t hash(char *key, uint32_t len)
{
    uint32_t   hash, i;
    for (hash=0, i=0; i<len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
} 

// SamplerRendererTask Definitions
void SamplerRendererTask::Run() {
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _SamplerRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }

    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];

    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
            float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

            // Evaluate radiance along camera ray
            PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
            if (visualizeObjectIds) {
                if (rayWeight > 0.f && scene->Intersect(rays[i], &isects[i])) {
                    // random shading based on shape id...
                    uint32_t ids[2] = { isects[i].shapeId, isects[i].primitiveId };
                    uint32_t h = hash((char *)ids, sizeof(ids));
                    float rgb[3] = { (h & 0xff), (h >> 8) & 0xff, (h >> 16) & 0xff };
                    Ls[i] = Spectrum::FromRGB(rgb);
                    Ls[i] /= 255.f;
                }
                else
                    Ls[i] = 0.f;
            }
            else {
            if (rayWeight > 0.f)
                Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
                                                 arena, &isects[i], &Ts[i]);
            else {
                Ls[i] = 0.f;
                Ts[i] = 1.f;
            }

            // Issue warning if unexpected radiance value returned
            if (Ls[i].HasNaNs()) {
                Error("Not-a-number radiance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            else if (Ls[i].y() < -1e-5) {
                Error("Negative luminance value, %f, returned"
                      "for image sample.  Setting to black.", Ls[i].y());
                Ls[i] = Spectrum(0.f);
            }
            else if (isinf(Ls[i].y())) {
                Error("Infinite luminance value returned"
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            }
            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
        }

        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, isects, sampleCount))
        {
            for (int i = 0; i < sampleCount; ++i)
            {
                PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
                camera->film->AddSample(samples[i], Ls[i]);
                PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
            }
        }

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }

    // Clean up after _SamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart,
        sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}



// SamplerRenderer Method Definitions
SamplerRenderer::SamplerRenderer(Sampler *s, Camera *c,
                                 SurfaceIntegrator *si, VolumeIntegrator *vi,
                                 bool visIds) {
    sampler = s;
    camera = c;
    surfaceIntegrator = si;
    volumeIntegrator = vi;
    visualizeObjectIds = visIds;
}


SamplerRenderer::~SamplerRenderer() {
    delete sampler;
    delete camera;
    delete surfaceIntegrator;
    delete volumeIntegrator;
}


void SamplerRenderer::Render(const Scene *scene) {
	    
	PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();

    surfaceIntegrator->Preprocess(scene, camera, this);
    volumeIntegrator->Preprocess(scene, camera, this);
    PBRT_FINISHED_PREPROCESSING();


    PBRT_STARTED_RENDERING();

if (!PhotonImage) {
    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, surfaceIntegrator,
                                volumeIntegrator, scene);

    // Create and launch _SamplerRendererTask_s for rendering image

    // Compute number of _SamplerRendererTask_s to create for rendering
    int nPixels = camera->film->xResolution * camera->film->yResolution;
    int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
    nTasks = RoundUpPow2(nTasks);
    ProgressReporter reporter(nTasks, "Rendering");
    vector<Task *> renderTasks;



    for (int i = 0; i < nTasks; ++i)
        renderTasks.push_back(new SamplerRendererTask(scene, this, camera,
                                                      reporter, sampler, sample, 
                                                      visualizeObjectIds, 
                                                      nTasks-1-i, nTasks));
    EnqueueTasks(renderTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < renderTasks.size(); ++i)
        delete renderTasks[i];
    reporter.Done();
    PBRT_FINISHED_RENDERING();
    // Clean up after rendering and store final image
    delete sample;
    camera->film->WriteImage();
}
}


Spectrum SamplerRenderer::Li(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect))
        Li = surfaceIntegrator->Li(scene, this, ray, *isect, sample,
                                   rng, arena);
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           Li += scene->lights[i]->Le(ray);
    }



Spectrum Lvi=0;
	if (PhotonImage)	//AJOUT POUR FAIRE DE L'ABSORPTION VOLUME
{

//    Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
  //                                      T, arena);
	
Lvi=0;
	*T=1;

}
else { 
Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
                                        T, arena);

}
    return *T * Li + Lvi;
}


Spectrum SamplerRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}


