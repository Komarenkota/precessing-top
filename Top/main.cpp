#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "top1.cpp" 
#include <fstream>

void drawTop(const Top& top) {
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    GLfloat rotation[16] = {
        static_cast<GLfloat>(top.R[0][0]), static_cast<GLfloat>(top.R[1][0]), static_cast<GLfloat>(top.R[2][0]), 0,
        static_cast<GLfloat>(top.R[0][1]), static_cast<GLfloat>(top.R[1][1]), static_cast<GLfloat>(top.R[2][1]), 0,
        static_cast<GLfloat>(top.R[0][2]), static_cast<GLfloat>(top.R[1][2]), static_cast<GLfloat>(top.R[2][2]), 0,
        0, 0, 0, 1
    };
    
    glPushMatrix();
    
    glTranslatef(top.x[0] , top.x[2], top.x[1]);
    glRotatef(-90.0f, 1.0f, 0.0f, 0.0f); 
    glMultMatrixf(rotation);

    glBegin(GL_LINES);
    glColor3f(1,0,0); glVertex3f(0,0,0); glVertex3f(0.3,0,0); // X`
    glColor3f(0,1,0); glVertex3f(0,0,0); glVertex3f(0,0.3,0); // Y`
    glColor3f(0,0,1); glVertex3f(0,0,0); glVertex3f(0,0,0.3); // Z`
    glEnd();

    GLUquadric* quadric = gluNewQuadric();
    glColor3f(1.0f, 0.0f, 0.0f);
    gluDisk(quadric, 0.0, 0.05, 32, 1);
    
    glTranslatef(top.r0[0], top.r0[1],top.r0[2]);
    glColor3f(0.0f, 1.0f, 0.0f);
    gluCylinder(quadric, 0.0, 0.05, 0.1, 32, 32);


    gluDeleteQuadric(quadric);
    glPopMatrix();
}


int main(int argc, char* argv[]) {
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("Spinning Top Simulation",
                                        SDL_WINDOWPOS_CENTERED,
                                        SDL_WINDOWPOS_CENTERED,
                                        800, 600,
                                        SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    
    SDL_GLContext context = SDL_GL_CreateContext(window);
    
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    
    glMatrixMode(GL_PROJECTION);
    gluPerspective(45.0, 800.0/600.0, 0.1, 100.0);
    glMatrixMode(GL_MODELVIEW);
    
    double T = 1000.0;
    double dt = 0.0001; 
    Top top = {0.2, {{0.5, 0, 0}, {0, 0.5, 0}, {0, 0, 0.05}}, {{2, 0, 0}, {0, 2, 0}, {0, 0, 20}}, {0.0, 0.0, -0.1}};

    initState(&top);
    double x0[13], x_fin[13];
    computeTorqueAndForces(&top);
    StateToArray(&top, x0);
    
    SDL_Event event;
    bool running = true;
    double t = 0;
    std::ofstream outx("Lx.txt");
    std::ofstream outy("Ly.txt");
    std::ofstream outxy("Lxy.txt");
    while (running && t < T) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
        }

        rk4(&top,x0, dt, x_fin);
        computeTorqueAndForces(&top);

        


        printf("V: %+.15le %+.15le %+.15le\n", top.v[0], top.v[1], top.v[2]);
        printf("F: %+.15le %+.15le %+.15le\n", top.F[0], top.F[1], top.F[2]);        
        printf("L: %+.15le %+.15le %+.15le\n", top.L[0], top.L[1], top.L[2]);
        printf("Xz: %+.15le \n", top.x[2]);
        // Вывод угла наклона в градусах
        double axis_world[3];
        multiplyMatrixByVector(top.R, top.r0, axis_world);
        double tilt_angle = computeTiltAngle(axis_world);
        std::cout << "t = " << t << ", Угол наклона = " << 180.0 - tilt_angle * 180.0 / M_PI << "°\n";
        checkEnergy(&top);

        outx << t << " "<<top.L[0] << '\n';

        
        outy << t <<" " <<top.L[1] << '\n';

        outxy << top.L[0] <<" " <<top.L[1] << '\n';

        for (int i = 0; i < 13; i++) x0[i] = x_fin[i];
        t += dt;
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glLoadIdentity();
        gluLookAt(1.0, 1.0, 1.0,
                 0.0, 0.0, 0.0,   
                 0.0, 0.1, 0.0);  

        
        
        drawTop(top);
        
        glBegin(GL_LINES);
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(1.0f, 0.0f, 0.0f);
        
        glColor3f(0.0f, 1.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, 1.0f, 0.0f);
        
        glColor3f(0.0f, 0.0f, 1.0f);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, 1.0f);
        glEnd();
        
        SDL_GL_SwapWindow(window);
        SDL_Delay(16); 
    }
    
    SDL_GL_DeleteContext(context);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}