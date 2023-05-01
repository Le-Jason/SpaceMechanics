import sys
sys.path.append('C:\\Users\\spong\\Documents\\Code\\projects\\manimrepo\\manim')
from manimlib import *
import numpy as np
sys.path.append('C:/Users/spong/Documents/Code/projects/spacemechanics/tools/data')
sys.path.append('C:/Users/spong/Documents/Code/projects/spacemechanics/tools')
from planetaryData import earth
from spaceCraftState import spaceCraftState
from spaceCraftState import perturbations
import math as m

class Orbits(Scene):
    def construct(self):


        r_apo = 7
        r_peri = 3
        a = (r_apo + r_peri) / 2
        e = (r_apo - r_peri)/(r_apo + r_peri)

        state = [a*2000,e,0,0,0,0]
        perts=perturbations()
        perts['J2'] + True
        Orbit1 = spaceCraftState('Earth',kepler=True,stateVec=state,perturbation=perts)
        t0 = 0
        y0 = [0,0,0,0,0,0]
        dt = 100
        tf = 100*5000
        tSol,ySol = Orbit1.propagateOrbit(t0,y0,dt,tf)

        

        intro_words = Text("""
            Fundamentals of Orbital Mechanics: 

                          Orbits
        """).scale(1)
        intro_words.to_edge(UP)
        orbit = Ellipse(width=10.0, height=4.0).shift(DOWN*1)



        circle = Circle()
        orbit.set_color(YELLOW)
        img = ImageMobject("Earth_Western_Hemisphere_transparent_background.png")
        img.set_height(circle.get_height())
        img.set_width(circle.get_width())
        img.shift(DOWN*1).shift(RIGHT*2)
        # print(orbit.get_anchors())
        maxValX = max(-1*ySol[:,0]) 
        maxValY = max(ySol[:,1]) 

        img2 = ImageMobject("Transiting_Exoplanet_Survey_Satellite_artist_concept_(transparent_background).png")
        img2 = Circle().scale(0.1)
        img3 = Circle().scale(0.1)
        img3.shift(LEFT*0).shift(UP*1.0)
        # self.add(img3)
        img2.p = ySol[0,0:3]
        self.t=0
        self.step=0.001
        self.counter=0
        def update(mob):
            mob.move_to([img.get_center()[0] + ySol[self.counter,0]*7/maxValX,img.get_center()[1] + ySol[self.counter,1]*2/maxValY,img.get_center()[2] + ySol[self.counter,2]])
            img2.p = [img.get_center()[0] + ySol[self.counter,0]*7/maxValX,img.get_center()[1] + ySol[self.counter,1]*2/maxValY,img.get_center()[2] + ySol[self.counter,2]]
            self.counter+=1
        self.add(orbit,img,img2)
        img2.add_updater(update)
        
        path=Line(img2.get_center(),img.get_center()).set_stroke(width=0.05)
        def path_update(mob):
            line=Line(img2.get_center(),img.get_center()).set_stroke(width=0.05)
            # path.append_vectorized_mobject(line)
            mob.become(line)

        def update_label(mob):
            pos = path.point_from_proportion(0.25)
            mob.move_to(pos + [0,0.5,0])
        

        path.add_updater(path_update)
        self.add(path)
        text2 = Text("r").scale(1)
        text2.add_updater(update_label)
        
        refAngle = path.point_from_proportion(0.25)/np.linalg.norm(path.point_from_proportion(0.25))
        scaleRelLine = 1.5
        refLine = path.point_from_proportion(1.0)- (scaleRelLine*path.get_unit_vector())
        refLineUnit = path.get_unit_vector()
        cheese = path.get_length()
        refLength = np.linalg.norm(path.point_from_proportion(0.6))
        self.oldAngle = 0
        self.CNT = 0
        def update_arc(mob):
            
            referenceLine = img2.get_center() - img.get_center()
            length = np.linalg.norm(referenceLine)
            coords = cheese/path.get_length()
            line2 = path.point_from_proportion(1.0) - (scaleRelLine*path.get_unit_vector())
            
            dot_pro = np.dot(refLineUnit, path.get_unit_vector())
            angle = m.degrees(m.acos(dot_pro))
            
            if angle < self.oldAngle:
                self.CNT = 1
            self.oldAngle = angle

            if self.CNT == 1:
                angle = (180 - angle + 180)
                self.CNT = 0
            if angle >= 320:
                angle = 320
            arc2 = ArcBetweenPoints(refLine,line2,angle*3.1415/180, stroke_color=WHITE,radius=1)
            mob.become(arc2)
        arc = ArcBetweenPoints(path.get_unit_vector(),path.get_unit_vector(), stroke_color=WHITE,radius=1)
        arc.add_updater(update_arc)
        def update_label2(mob):
            pos = path.point_from_proportion(1.0) - (0.75*scaleRelLine*path.get_unit_vector())
            mob.move_to(pos)


        def update_path2(mob):
            theta = -np.pi/2
            R = [m.cos(theta), -m.sin(theta), 0,
                 m.sin(theta), m.cos(theta), 0,
                 0, 0, 1]
            R = np.array(R)
            R = R.reshape((3,3))
            
            rotLine = np.matmul(R,path.get_unit_vector())
            print(R)
            line3=Line(path.get_start(),path.get_unit_vector()).set_stroke(width=0.50)
            mob.become(line3)
        text3 = Text("v").scale(1)
        text3.add_updater(update_label2)
        path2 = Line(img2.get_center(),img2.get_center()).set_stroke(width=0.50)
        path2.add_updater(update_path2)
    
        self.add(path2)
        self.add(text3)
        
        self.add(arc)
        self.add(text2)
        self.play(Write(intro_words))
        self.wait(60)
