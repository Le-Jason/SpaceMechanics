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

        # Setting up the Orbit
        r_apo = 7 # Change the Width of the orbit
        r_peri = 3 # Change the Height of the orbit
        orbitSpeed = 20000 # Change the Speed of the orbit

        # Create the Orbit
        a = (r_apo + r_peri) / 2
        e = (r_apo - r_peri)/(r_apo + r_peri)
        state = [a*orbitSpeed,e,0,0,0,0] 
        perts=perturbations()
        perts['J2'] + True
        Orbit1 = spaceCraftState('Earth',kepler=True,stateVec=state,perturbation=perts)
        t0 = 0
        y0 = [0,0,0,0,0,0]
        dt = 100
        tf = 100*5000
        tSol,ySol = Orbit1.propagateOrbit(t0,y0,dt,tf)
        

        # Setup Orbit's horizonal and vertical lines
        periCounter = ySol[0,0]
        periPos = ySol[0,0:3]
        apogeeCounter = ySol[0,0]
        apogeePos = ySol[0,0:3]
        yAxisHighCounter = ySol[0,1]
        yAxisHighPos = ySol[0,0:3]
        yAxisLowCounter = ySol[0,1]
        yAxisLowPos = ySol[0,0:3]
        for i in range(len(ySol)):
            x = ySol[i,0]
            y = ySol[i,1]
            z = ySol[i,2]
            if periCounter > x:
                periCounter = x
                periPos = ySol[i,0:3]
            if apogeeCounter < x:
                apogeeCounter = x
                apogeePos = ySol[i,0:3]
            if yAxisLowCounter < y:
                yAxisLowCounter = y
                yAxisLowPos = ySol[i,0:3]
            if yAxisHighCounter > y:
                yAxisHighCounter = y
                yAxisHighPos = ySol[i,0:3]
        maxValX = max(-1*ySol[:,0]) 
        maxValY = max(ySol[:,1]) 

        # Scene SetUp
        intro_words = Text("""
            Fundamentals of Orbital Mechanics: 

                          Orbits
        """).scale(1)
        intro_words.to_edge(UP)

        # Planet
        orbit = Ellipse(width=10.0, height=4.0).shift(DOWN*1)
        circle = Circle()
        orbit.set_color(YELLOW)
        img = ImageMobject("Earth_Western_Hemisphere_transparent_background.png")
        img.set_height(circle.get_height())
        img.set_width(circle.get_width())
        img.shift(DOWN*1).shift(RIGHT*2)
        
        # Orbit
        
        yAxixLine = [yAxisLowPos[0]*7/maxValX,yAxisLowPos[1]*15/5/maxValX,yAxisLowPos[2]*7/maxValX]
        yAxixLine2 = [yAxisHighPos[0]*7/maxValX,yAxisHighPos[1]*15/5/maxValX,yAxisHighPos[2]*7/maxValX]
        orbitH=Line(img.get_center()[:] + apogeePos*7/maxValX,img.get_center()[:] + periPos*7/maxValX).set_stroke(width=0.05)
        orbitV=Line(img.get_center()[:] + yAxixLine,img.get_center()[:] +yAxixLine2).set_stroke(width=0.05)
        
        
        

        img2 = ImageMobject("surveysatellite.png")
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
            line3=Line(-rotLine+path.get_start(),rotLine+path.get_start()).set_stroke(width=0.50)
            mob.become(line3)
        path2 = Line(img2.get_center(),img2.get_center()).set_stroke(width=0.50)
        path2.add_updater(update_path2)

        
        def update_arc4(mob):
            linePath1 = [img.get_center()[0] + ySol[self.counter-1,0]*7/maxValX,img.get_center()[1] + ySol[self.counter-1,1]*2/maxValY,img.get_center()[2] + ySol[self.counter-1,2]]
            linePath2 = [img.get_center()[0] + ySol[self.counter+1,0]*7/maxValX,img.get_center()[1] + ySol[self.counter+1,1]*2/maxValY,img.get_center()[2] + ySol[self.counter+1,2]]
            intermedLine=Line(linePath2,linePath1).set_stroke(width=0.50)
            refLineRef = intermedLine.point_from_proportion(1.0)- (0.5*intermedLine.get_unit_vector())
            line2Ref = path2.point_from_proportion(1.0) - (0.5*path2.get_unit_vector())
            dot_pro = np.dot(intermedLine.get_unit_vector(), path2.get_unit_vector())
            angle = m.degrees(m.acos(dot_pro))
            
            cross_angle = np.cross(intermedLine.get_unit_vector(),path2.get_unit_vector())
            if cross_angle[2] > 0:
                arc3 = ArcBetweenPoints(line2Ref,refLineRef,angle*3.1415/180/2, stroke_color=WHITE,radius=1)
            else:
                arc3 = ArcBetweenPoints(refLineRef,line2Ref,angle*3.1415/180/2, stroke_color=WHITE,radius=1)
                
            mob.become(arc3)
        def update_line4(mob):
            linePath1 = [img.get_center()[0] + ySol[self.counter-1,0]*7/maxValX,img.get_center()[1] + ySol[self.counter-1,1]*2/maxValY,img.get_center()[2] + ySol[self.counter-1,2]]
            linePath2 = [img.get_center()[0] + ySol[self.counter+1,0]*7/maxValX,img.get_center()[1] + ySol[self.counter+1,1]*2/maxValY,img.get_center()[2] + ySol[self.counter+1,2]]
            intermedLine=Line(linePath1,linePath2).set_stroke(width=0.50)
            refLineRef = intermedLine.point_from_proportion(1.0)- (0.5*intermedLine.get_unit_vector())
            line2Ref = path2.point_from_proportion(1.0) - (0.5*path2.get_unit_vector())
            dot_pro = np.dot(intermedLine.get_unit_vector(), path2.get_unit_vector())
            angle = m.degrees(m.acos(dot_pro))
            cross_angle = np.cross(intermedLine.get_unit_vector(),path2.get_unit_vector())
            theta1 = -np.pi/2
            R1 = [m.cos(theta1), -m.sin(theta1), 0,
                 m.sin(theta1), m.cos(theta1), 0,
                 0, 0, 1]
            R1 = np.array(R1)
            R1 = R1.reshape((3,3))
            if cross_angle[2] > 0:
                theta2 = -m.acos(dot_pro)
            else:
                theta2 = m.acos(dot_pro)
            
            R2 = [m.cos(theta2), -m.sin(theta2), 0,
                 m.sin(theta2), m.cos(theta2), 0,
                 0, 0, 1]
            R2 = np.array(R2)
            R2 = R2.reshape((3,3))
            
            rotLine = np.matmul(R1,path.get_unit_vector())

            rotLine = np.matmul(R2,rotLine)
            line4=Line(-0+path.get_start(),2*rotLine+path.get_start()).set_stroke(width=4)
            mob.become(line4)
        path4 = Line([img.get_center()[0] + ySol[self.counter-1,0]*7/maxValX,img.get_center()[1] + ySol[self.counter-1,1]*2/maxValY,img.get_center()[2] + ySol[self.counter-1,2]],[img.get_center()[0] + ySol[self.counter+1,0]*7/maxValX,img.get_center()[1] + ySol[self.counter+1,1]*2/maxValY,img.get_center()[2] + ySol[self.counter+1,2]]).set_stroke(width=0.50)
        path4.add_updater(update_line4)
        
        arc4 = ArcBetweenPoints(path.get_unit_vector(),path.get_unit_vector(), stroke_color=WHITE,radius=1)
        arc4.add_updater(update_arc4)
        text3 = Text("v").scale(1)
        text3.add_updater(update_label2)
        path4.add_updater(update_line4)
    
        self.add(path2)
        # self.add(text3)
        self.add(path4)
        self.add(arc)
        self.add(arc4)
        # self.add(text2)
        self.add(orbitH)
        self.add(orbitV)
        self.play(Write(intro_words))
        self.wait(60)