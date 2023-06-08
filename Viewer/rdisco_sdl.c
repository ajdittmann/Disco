#include <stdio.h>
#include <math.h>
#include <SDL2/SDL.h>
#include "rdisco.h"

const int INIT_SCREEN_WIDTH = 600;
const int INIT_SCREEN_HEIGHT = 400;

void gaussian_func(double xyz[3], double *prim, int Nq, float rgba[4])
{
    int N = 6;
    float rhos[] = {1.0e-5, 1.0e-4, 1.0e-3, 2.0e-2, 3.0e-2, 1.0};
    float alphas[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    float widthFac[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

    float cols[] = {0.0, 0.0, 1.0,
                    0.0, 1.0, 1.0,
                    0.0, 1.0, 0.0,
                    1.0, 1.0, 0.0,
                    1.0, 0.0, 0.0,
                    1.0, 1.0, 1.0};

    float r = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    float z = xyz[2];

    double temp = prim[1] / prim[0];

    rgba[0] = 1.0*temp * prim[0]; //1.0;
    rgba[1] = 40*temp*temp * prim[0]; //1.0;
    rgba[2] = 5e4*temp*temp*temp*temp * prim[0]; //1.0;
    rgba[3] = 50.0*prim[0];

    return;
    rgba[0] = 0.0;
    rgba[1] = 0.0;
    rgba[2] = 0.0;
    rgba[3] = 0.0;

    int i;
    for(i=2; i<N; i++)
    {
        float x = ((float)(prim[0]) - rhos[i]) / (widthFac[i] * rhos[i]);
        if(fabs(x) > 5)
            continue;

        float f = exp(-0.5*x*x);

        rgba[0] += f*cols[3*i + 0];
        rgba[1] += f*cols[3*i + 1];
        rgba[2] += f*cols[3*i + 2];
        rgba[3] += f*alphas[i]*1.0e-3;
    }
}

void normalize_image(float *im, int Nx, int Ny)
{
    printf("Normalizing image\n");
    int i;
    float max_channel = 0;
    for(i=0; i<3*Nx*Ny; i++)
    {
        if(im[i] > max_channel)
            max_channel = im[i];
    }

    float im_norm = 1.0 / max_channel;

    im_norm = 1.4e4;
    printf("im_norm: %.3g\n", im_norm);

    for(i=0; i<3*Nx*Ny; i++)
        im[i] *= im_norm;
}

void copyImageToSurface(float *im, int Nx, int Ny, SDL_Surface *surface)
{
    printf("Copying image to SDL Surface\n");
    SDL_LockSurface(surface);

    int i, j;
    
    for(i=0; i<Nx; i++)
        for(j=0; j<Ny; j++)
        {
            int img_idx = 3 * ((Ny-j-1) * Nx + i);
            unsigned int r = 255.99 * im[img_idx  ];
            unsigned int g = 255.99 * im[img_idx + 1];
            unsigned int b = 255.99 * im[img_idx + 2];

            Uint32 *pix = (Uint32 *) ((Uint8 *)surface->pixels
                            + surface->pitch * j
                            + surface->format->BytesPerPixel * i);
            *pix = SDL_MapRGBA(surface->format, r, g, b, 255);
        }
    
    SDL_UnlockSurface(surface);
}

int main(int argc, char *argv[])
{
    SDL_Window *window = NULL;

    if(SDL_Init(SDL_INIT_VIDEO) < 0 )
    {
        printf("SDL could not initialize! SDL Error: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 0;
    }

    window = SDL_CreateWindow("Render Disco",
                             SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                             INIT_SCREEN_WIDTH, INIT_SCREEN_HEIGHT,
                             SDL_WINDOW_SHOWN);

    if (window == NULL)
    {
        printf("Window could not be created! SDL_Error %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 0;
    }

    int width, height;

    SDL_GetWindowSize(window, &width, &height);

    SDL_Surface *windowSurface;
    windowSurface = SDL_GetWindowSurface(window);
    if(windowSurface == NULL)
    {
        printf("Window Surface could not be created! SDL_Error %s\n",
                SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 0;
    }

    SDL_Surface *imgSurface;
    imgSurface = SDL_CreateRGBSurface(0, width, height,
                                   32, 0xFF000000, 0x00FF0000, 0x0000FF00,
                                   0x000000FF);
    if(imgSurface == NULL)
    {
        printf("Image Surface could not be created! SDL_Error %s\n",
                SDL_GetError());
        SDL_FreeSurface(imgSurface);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 0;
    }

    int file_idx = 0;
    char **filenames = argv + 1;
    int n_files = argc - 1;
    
    printf("Loading file: %s\n", filenames[file_idx]);
    struct Grid g = grid_new(filenames[file_idx]);
    grid_print(&g);

    //double cam_xyz[3] = {0.0, 0.0, 20.0};
    //double cam_xyz[3] = {0.0, -17.0, 10.0};
    //double cam_nx[3] = {1.0, 0.0, 0.0};
    double cam_xyz[3] = {17.0, 0.0, 10.0};
    double cam_nx[3] = {0.0, 1.0, 0.0};

    cam_xyz[0] /= 2;
    cam_xyz[1] /= 2;
    cam_xyz[2] /= 2;
    
    double origin[3] = {0.0, 0.0, 0.0};
    double cam_nz[3];
    look_at(origin, cam_xyz, cam_nz);

    struct Ray *rays;

    //rays = ray_generate_ortho(width, height, cam_xyz, cam_nz, cam_nx, 8.0);
    rays = ray_generate_persp(width, height, cam_xyz, cam_nz, cam_nx, 60.0);

    struct IntersectionList *ray_hits = (struct IntersectionList *) malloc(
            width * height * sizeof(struct IntersectionList));

    int i;
    printf("Intersecting...\n");
    
    for(i=0; i<width * height; i++)
    {
        ray_hits[i] = intersect(&g, rays+i);
        /*
        printf("Ray %d: %d hits (%g -> %g)\n",
                i, ray_hits[i].len,
                ray_hits[i].buf[0].sin,
                ray_hits[i].buf[ray_hits[i].len-1].sout);
        */
    }
        

    printf("Intersection Done.\n");

    float *image = (float *)malloc(width*height*3*sizeof(float));

    render(&g, rays, ray_hits, image, width, height, gaussian_func);
    
    printf("Render Done.\n");
    
    normalize_image(image, width, height);
    
    printf("Scaling Done.\n");

    copyImageToSurface(image, width, height, imgSurface);

    int err = SDL_BlitSurface(imgSurface, NULL, windowSurface, NULL);
    if(err)
    {
        printf("Image could not be blitted! SDL_Error %s\n",
                SDL_GetError());
        SDL_FreeSurface(imgSurface);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 0;
    }

    SDL_UpdateWindowSurface(window);

    int quit = 0;

    while(!quit)
    {
        SDL_Event e;
        while(SDL_PollEvent(&e))
        {
            if(e.type == SDL_QUIT)
                quit = 1;
            else if(e.type == SDL_KEYDOWN)
            {
                switch(e.key.keysym.sym)
                {
                    case SDLK_ESCAPE:
                        quit = 1;
                        break;
                    case SDLK_COMMA:
                        file_idx--;
                        if(file_idx < 0)
                            file_idx = n_files - 1;
                        grid_load(&g, filenames[file_idx]);
                        refreshIntersections(&g, rays, ray_hits, width*height);
                        render(&g, rays, ray_hits, image, width, height,
                                gaussian_func);
                        normalize_image(image, width, height);
                        copyImageToSurface(image, width, height, imgSurface);
                        SDL_BlitSurface(imgSurface, NULL, windowSurface, NULL);
                        SDL_UpdateWindowSurface(window);
                        break;
                    case SDLK_PERIOD:
                        file_idx++;
                        if(file_idx >= n_files)
                            file_idx = 0;
                        grid_load(&g, filenames[file_idx]);
                        refreshIntersections(&g, rays, ray_hits, width*height);
                        render(&g, rays, ray_hits, image, width, height,
                                gaussian_func);
                        normalize_image(image, width, height);
                        copyImageToSurface(image, width, height, imgSurface);
                        SDL_BlitSurface(imgSurface, NULL, windowSurface, NULL);
                        SDL_UpdateWindowSurface(window);
                        break;
                    case SDLK_m:
                        for(file_idx=0; file_idx<n_files; file_idx++)
                        {
                            grid_load(&g, filenames[file_idx]);
                            refreshIntersections(&g, rays, ray_hits,
                                                  width*height);
                            render(&g, rays, ray_hits, image, width, height,
                                    gaussian_func);
                            normalize_image(image, width, height);
                            copyImageToSurface(image, width, height,
                                                imgSurface);
                            SDL_BlitSurface(imgSurface, NULL,
                                            windowSurface, NULL);
                            SDL_UpdateWindowSurface(window);

                            char imgname[256];
                            sprintf(imgname, "%s.bmp", filenames[file_idx]);
                            SDL_SaveBMP(imgSurface, imgname);
                        }
                        break;
                    case SDLK_s:
                        
                        printf("Saving\n");
                        char imgname[256];
                        sprintf(imgname, "%s.bmp", filenames[file_idx]);
                        SDL_SaveBMP(imgSurface, imgname);
                        printf("Saved to %s\n", imgname);
                        
                        break;
                }
            }
        }
    }
    
    SDL_FreeSurface(imgSurface);

    free(image);

    for(i=0; i<width*height; i++)
        intersectList_free(ray_hits + i);
    free(ray_hits);

    free(rays);
    grid_free(&g);

    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
