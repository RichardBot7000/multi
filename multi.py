from big_ol_pile_of_manim_imports import *
import numpy as np
import itertools as it
from copy import deepcopy
import sys
from manimlib.constants import *
from manimlib.scene.scene import Scene
from manimlib.mobject.geometry import Polygon
from manimlib.once_useful_constructs.region import  region_from_polygon_vertices, region_from_line_boundary


const = [
	0.1, #constant in the exponent
]

time = 0

class TestSurface(ThreeDScene):
	def construct(self):
		self.set_camera_orientation(phi=60*DEGREES, theta = 30*DEGREES)
		surface = ParametricSurface(
			lambda u, v: np.array([
				u,
				v,
				3*np.sin(u)*np.sin(v)/u/v
			]),v_min=-5,v_max=5,u_min=-5,u_max=5,
			resolution=(25, 25)).fade(0.5)
		surface.color_using_background_image("/Users/richardwu/Documents/Manim/manim_3feb/images/TG.jpeg")
		self.play(Write(surface))
		self.move_camera(theta=120*DEGREES,run_time=3)
		# self.move_camera(phi=-15*DEGREES, run_time = 1)



class FirstEquation(GraphScene):
	CONFIG = {
		"rate_func": linear,
		"y_axis_label": "$T$",
	}
	def construct(self):
		self.setup_axes()
		dot = Dot().scale(0.00001)
		def func(x):
			return 2*np.sin(2*PI*x/3)*(np.e**(-const[0]*((2*PI/3)**2)*dot.get_center()[0])) + 2*np.sin(2*PI*x/5)*(np.e**(-const[0]*((2*PI/5)**2)*dot.get_center()[0])) + np.sin(2*PI*x/7)*(np.e**(-const[0]*((2*PI/7)**2)*dot.get_center()[0]))+4
		input_tracker_1 = ValueTracker(6)
		graph = self.get_graph(func,x_min=0,x_max=9.191)
		graph.color_using_background_image("/Users/richardwu/Documents/Manim/manim_3feb/images/TG.jpeg")
		def get_x_value(input_tracker):
			return input_tracker.get_value()
		def get_y_value(input_tracker):
			return graph.underlying_function(get_x_value(input_tracker))
		def get_graph_point(input_tracker):
			return self.coords_to_point(get_x_value(input_tracker), get_y_value(input_tracker))
		eqn = [
			TexMobject("T(x", ")", "=", "2\\sin(\\frac{2}{3}\\pi x) + 2\\sin(\\frac{2}{5}\\pi x) + \\sin(\\frac{2}{7}\\pi x) + 4").move_to([0.5, 2.8, 0]).scale(0.7),
			TexMobject("T(x", ", 0", ")", "=", "2\\sin(\\frac{2}{3}\\pi x) + 2\\sin(\\frac{2}{5}\\pi x) + \\sin(\\frac{2}{7}\\pi x) + 4").move_to([0.5, 2.8, 0]).scale(0.7)
		]
		colorarray = [
			[WHITE, WHITE, WHITE, WHITE],
			[WHITE, WHITE, WHITE, WHITE, WHITE],
		]
		for i in range(len(eqn)):
			for j in range(len(eqn[i])):
				eqn[i][j].set_color(colorarray[i][j])
		self.play(ShowCreation(graph), run_time = 2)
		self.play(Write(eqn[0]))
		self.wait()
		changes = [
			[[(0, 1, 2, 3), (0, 2, 3, 4)]],
			[[(0, 1, 2, 3, 4, 1), (0, 2, 1, 4, 5, 3)]]
		]
		for pre_ind, post_ind in changes[0]:
			self.play(*[Transform(eqn[0][i], eqn[1][j]) for i,j in zip(pre_ind, post_ind)], FadeIn(eqn[1][1]))
		self.wait(2)
		text = TexMobject("time = ").move_to([5, 0 ,0])
		timer = DecimalNumber(0, num_decimal_places = 2).next_to(text, RIGHT, buff = 0.2)
		self.play(FadeIn(timer), FadeIn(text))
		group = VGroup(timer, graph)
		def update_timer(group):
			timer = DecimalNumber(dot.get_center()[0], num_decimal_places = 2).next_to(text, RIGHT, buff = 0.2)
			graph = self.get_graph(func, x_min = 0, x_max = 9.191)
			graph.color_using_background_image("/Users/richardwu/Documents/Manim/manim_3feb/images/TG.jpeg")
			new_group = VGroup(timer, graph)
			group.become(new_group)
			return group
		self.wait()
		self.play(dot.shift, 10*RIGHT, 
			UpdateFromFunc(group, update_timer), 
			rate_func = linear,
			run_time = 10)
		self.wait()


class TangentLine(GraphScene):
	CONFIG = {
		"y_axis_label": "$T$",
	}
	def construct(self):
		self.setup_axes()
		dot = Dot().scale(0.00001)
		def func(x):
			return 2*np.sin(2*PI*x/3)*(np.e**(-const[0]*((2*PI/3)**2)*dot.get_center()[0])) + 2*np.sin(2*PI*x/5)*(np.e**(-const[0]*((2*PI/5)**2)*dot.get_center()[0])) + np.sin(2*PI*x/7)*(np.e**(-const[0]*((2*PI/7)**2)*dot.get_center()[0]))+4
		def deriv(x):
			return 2*np.cos(2*PI*x/3)*(np.e**(-const[0]*((2*PI/3)**2)*dot.get_center()[0]))*(2*PI/3) + 2*np.cos(2*PI*x/5)*(np.e**(-const[0]*((2*PI/5)**2)*dot.get_center()[0]))*(2*PI/5) + np.cos(2*PI*x/7)*(np.e**(-const[0]*((2*PI/7)**2)*dot.get_center()[0]))*(2*PI/7)
		def angle_from_deriv(x):
			return np.arctan(deriv(x))
		dot2 = Dot().scale(0.0000001).move_to([0, 0, 0])
		x_val = dot2.get_center()[0]
		input_tracker_1 = ValueTracker(x_val)
		graph = self.get_graph(func,x_min=0,x_max=9.191)
		graph.color_using_background_image("/Users/richardwu/Documents/Manim/manim_3feb/images/TG.jpeg")
		def get_x_value(input_tracker):
			return input_tracker.get_value()
		def get_y_value(input_tracker):
			return graph.underlying_function(get_x_value(input_tracker))
		def get_graph_point(input_tracker):
			return self.coords_to_point(get_x_value(input_tracker), get_y_value(input_tracker))
		bot = Dot().move_to(get_graph_point(input_tracker_1))
		line = Line([bot.get_center()[0]-np.cos(angle_from_deriv(x_val)), bot.get_center()[1]-np.sin(angle_from_deriv(x_val)), 0], [bot.get_center()[0]+np.cos(angle_from_deriv(x_val)), bot.get_center()[1]+np.sin(angle_from_deriv(x_val)), 0])
		self.play(ShowCreation(graph), run_time = 2)
		text = TexMobject("time = ").move_to([4, 2 ,0])
		partial = TexMobject("\\frac{\\partial T}{\\partial x} = ").move_to([0, 3, 0])
		timer = DecimalNumber(0, num_decimal_places = 2).next_to(text, RIGHT, buff = 0.2)
		partial_deriv = DecimalNumber(deriv(x_val), num_decimal_places = 2).next_to(partial, RIGHT, buff = 0.2)
		self.play(FadeIn(timer), FadeIn(text))
		group = VGroup(timer, graph, bot, line, input_tracker_1, partial_deriv)
		def update_timer(group):
			timer = DecimalNumber(dot.get_center()[0], num_decimal_places = 2).next_to(text, RIGHT, buff = 0.2)
			graph = self.get_graph(func, x_min = 0, x_max = 9.191)
			graph.color_using_background_image("/Users/richardwu/Documents/Manim/manim_3feb/images/TG.jpeg")
			x_val = dot2.get_center()[0]
			input_tracker_1 = ValueTracker(x_val)
			bot = Dot().move_to(get_graph_point(input_tracker_1))
			line = Line([bot.get_center()[0]-np.cos(angle_from_deriv(x_val)), bot.get_center()[1]-np.sin(angle_from_deriv(x_val)), 0], [bot.get_center()[0]+np.cos(angle_from_deriv(x_val)), bot.get_center()[1]+np.sin(angle_from_deriv(x_val)), 0])
			partial_deriv = DecimalNumber(deriv(x_val), num_decimal_places = 2).next_to(partial, RIGHT, buff = 0.2)
			new_group = VGroup(timer, graph, bot, line, input_tracker_1, partial_deriv)
			group.become(new_group)
			return group
		self.play(FadeIn(bot))
		self.play(GrowFromCenter(line), FadeIn(partial), FadeIn(partial_deriv))
		self.wait()
		self.play(dot2.shift, 9*RIGHT, UpdateFromFunc(group, update_timer), rate_func = smooth, run_time = 8)
		self.play(dot2.shift, 3*LEFT, UpdateFromFunc(group, update_timer), rate_func = smooth, run_time = 3)
		self.wait(2)
		self.play(dot.shift, 10*RIGHT, 
			UpdateFromFunc(group, update_timer), 
			rate_func = linear,
			run_time = 10)
		self.wait()


class SolveEquation(Scene):
	def construct(self):
		eqn = [
			TexMobject("T(x, 0)", "=", "Asin(b(x-c))"),
			TexMobject("{\\partial", "T", "\\over", "\\partial t}", "=", "\\alpha", "{\\partial^2", "T", "\\over", "\\partial x^2}"),
			TexMobject("T(x, t)", "=", "Asin(b(x-c))e^{-\\alpha b^2t}"),
			TexMobject("{\\partial", "T", "\\over", "\\partial t}", "=", "\\alpha", "{\\partial^2", "T", "\\over", "\\partial x^2}"),
			TexMobject("{\\partial", "\\over", "\\partial t}", "T", "=", "\\alpha", "{\\partial^2", "\\over", "\\partial x^2}", "T"),
			TexMobject("{\\partial", "\\over", "\\partial t}", "Asin(b(x-c))e^{-\\alpha b^2t}", "=", "\\alpha", "{\\partial^2", "\\over", "\\partial x^2}", "Asin(b(x-c))e^{-\\alpha b^2t}"),
			TexMobject("{\\partial", "\\over", "\\partial t}", "Asin(b(x-c))", "e^{-\\alpha b^2t}", "=", "\\alpha", "{\\partial^2", "\\over", "\\partial x^2}", "Asin(b(x-c))", "e^{-\\alpha b^2t}"),
			TexMobject("-\\alpha b^2", "Asin(b(x-c))", "e^{-\\alpha b^2t}", "=", "\\alpha", "{\\partial^2", "\\over", "\\partial x^2}", "Asin(b(x-c))", "e^{-\\alpha b^2t}"),
			TexMobject("-\\alpha b^2", "Asin(b(x-c))", "e^{-\\alpha b^2t}", "=", "\\alpha", "b", "{\\partial", "\\over", "\\partial x}", "Acos(b(x-c))", "e^{-\\alpha b^2t}"),
			TexMobject("-\\alpha b^2", "Asin(b(x-c))", "e^{-\\alpha b^2t}", "=", "-", "\\alpha", "b^2", "Asin(b(x-c))", "e^{-\\alpha b^2t}"),
		]
		changes = [
			[[(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), (0, 3, 1, 2, 4, 5, 6, 9, 7, 8)]],
			[[(0, 1, 2, 4, 5, 6, 7, 8), (0, 1, 2, 4, 5, 6, 7, 8)]],
			[[(3, 4, 5, 6, 7, 8, 9, 10, 11, 4), (1, 2, 3, 4, 5, 6, 7, 8, 9, 0)]],
			[[(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), (0, 1, 2, 3, 4, 6, 7, 8, 9, 10)]],
			[[(0, 1, 2, 3, 4, 5, 9, 10, 7), (0, 1, 2, 3, 5, 6, 7, 8, 4)]],
		]
		for i in range(3):
			eqn[i].scale(0.9)
		for i in range(3, 10):
			eqn[i].scale(0.8)
		self.play(Write(eqn[0]))
		arr = [RED, ORANGE, RED, RED, YELLOW, GREEN, BLUE, PURPLE, BLUE, BLUE]
		arr1 = [RED, RED, RED, ORANGE, YELLOW, GREEN, BLUE, BLUE, BLUE, PURPLE]
		arr2 = [RED, RED, RED, ORANGE, ORANGE, YELLOW, GREEN, BLUE, BLUE, BLUE, PURPLE, PURPLE]
		arr3 = [[ORANGE, ORANGE, ORANGE, YELLOW, GREEN, BLUE, BLUE, BLUE, PURPLE, PURPLE],
			[ORANGE, ORANGE, ORANGE, YELLOW, GREEN, PURPLE, BLUE, BLUE, BLUE, PURPLE, PURPLE],
			[ORANGE, ORANGE, ORANGE, YELLOW, PURPLE, GREEN, PURPLE, PURPLE, PURPLE]]
		for i in range(len(arr)):
			eqn[1][i].set_color(arr[i])
			eqn[3][i].set_color(arr[i])
		for i in range(len(arr1)):
			eqn[4][i].set_color(arr1[i])
			eqn[5][i].set_color(arr1[i])
		for i in range(len(arr2)):
			eqn[6][i].set_color(arr2[i])
		for i in range(len(arr3)):
			for j in range(len(arr3[i])):
				eqn[7+i][j].set_color(arr3[i][j])
		LOC = [0, 3, 0]
		self.play(
			eqn[0].move_to, LOC,
			run_time = 2)
		eqn[1].next_to(eqn[0], DOWN)
		self.play(Write(eqn[1]))
		eqn[2].next_to(eqn[1], DOWN)
		self.play(Write(eqn[2]))
		box = SurroundingRectangle(eqn[2], color = GREEN)
		self.play(ShowCreation(box), run_time = 2)
		eqn[3].next_to(eqn[2], DOWN)
		eqn[4].move_to(eqn[3].get_center())
		eqn[5].move_to(eqn[4].get_center())
		for i in range(6, 9):
			eqn[i+1].next_to(eqn[i], DOWN)
		temp = eqn[1].copy()
		self.play(Transform(temp, eqn[3]), run_time = 2)
		self.add(eqn[3])
		self.remove(temp)
		for pre_ind, post_ind in changes[0]:
			self.play(*[Transform(eqn[3][i], eqn[4][j]) for i,j in zip(pre_ind, post_ind)])
		self.remove(eqn[3])
		self.add(eqn[4])
		tmp1 = eqn[2][2].copy()
		tmp2 = eqn[2][2].copy()
		for pre_ind, post_ind in changes[1]:
			self.play(*[Transform(eqn[4][i], eqn[5][j]) for i,j in zip(pre_ind, post_ind)], 
				Transform(tmp1, eqn[5][3]),
				Transform(tmp2, eqn[5][9]),
				FadeOut(eqn[4][3]), FadeOut(eqn[4][9]))
		self.remove(eqn[4], eqn[5])
		self.add(eqn[6])
		self.wait()
		for pre_ind, post_ind in changes[2]:
			self.play(*[Transform(eqn[6][i].copy(), eqn[7][j]) for i,j in zip(pre_ind, post_ind)])
		for pre_ind, post_ind in changes[3]:
			self.play(*[Transform(eqn[7][i].copy(), eqn[8][j]) for i,j in zip(pre_ind, post_ind)], FadeIn(eqn[8][5]))
		for pre_ind, post_ind in changes[4]:
			self.play(*[Transform(eqn[8][i].copy(), eqn[9][j]) for i,j in zip(pre_ind, post_ind)])
		self.wait()
		mob = TexMobject("\\checkmark").set_color("#19e625").next_to(eqn[9], UR, buff = 0.2)
		self.play(FadeIn(mob))
		self.wait(2)
		self.play(WiggleOutThenIn(eqn[2]), WiggleOutThenIn(box))
		self.wait(2)


class Fourier(GraphScene):
	CONFIG = {
		"y_axis_label": "$T$",
	}
	def construct(self):
		self.setup_axes()
		def func(x):
			return 0.25*(x-5)**2+1
		#was too lazy to learn how to input multiple variables, so enjoy the copy paste
		# def fourier(x, n):
		# 	ans = 37/12
		# 	for i in range(1, n+1):
		# 		ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
		# 	return ans
		def fourier0(x):
			ans = 37/12
			for i in range(1, 1):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		def fourier1(x):
			ans = 37/12
			for i in range(1, 2):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		def fourier2(x):
			ans = 37/12
			for i in range(1, 3):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		def fourier3(x):
			ans = 37/12
			for i in range(1, 4):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		def fourier4(x):
			ans = 37/12
			for i in range(1, 5):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		def fourier5(x):
			ans = 37/12
			for i in range(1, 6):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		def fourier10(x):
			ans = 37/12
			for i in range(1, 11):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		def fourier20(x):
			ans = 37/12
			for i in range(1, 21):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		def fourier50(x):
			ans = 37/12
			for i in range(1, 51):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		def fourier100(x):
			ans = 37/12
			for i in range(1, 101):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		def fourier1000(x):
			ans = 37/12
			for i in range(1, 1001):
				ans += ((-1)**(i))*25/(PI**2)/(i**2)*np.cos(i*PI*(x-5)/5)
			return ans
		graph = self.get_graph(func,x_min=0,x_max=10).color_using_background_image("/Users/richardwu/Documents/Manim/manim_3feb/images/TG.jpeg").set_stroke(width = 7)
		self.play(ShowCreation(graph), run_time = 2)
		val = ValueTracker(10)
		def get_x_value(input_tracker):
			return input_tracker.get_value()
		def get_y_value(input_tracker):
			return graph.underlying_function(get_x_value(input_tracker))
		def get_graph_point(input_tracker):
			return self.coords_to_point(get_x_value(input_tracker), get_y_value(input_tracker))
		graphlabel = TexMobject("f(x) = \\frac{1}{4}(x-5)^2+1").scale(0.6).next_to(get_graph_point(val), UP, buff = 0.1).shift(0.5*RIGHT)
		shiftedgraphlabel = TexMobject("s(x) = \\frac{1}{4}x^2+1").scale(0.6).next_to(get_graph_point(val), UP, buff = 0.1).shift(0.5*RIGHT)
		shiftedshiftedgraphlabel = graphlabel.copy()
		self.wait()
		self.play(FadeIn(graphlabel))
		LOC = [0, 3, 0]
		eqn = [
			TexMobject("s_N(x)", "=", "?").move_to(LOC),
			TexMobject("s_N(x)", "=", "\\frac{a_0}{2}", "+", "\\sum_{n=1}^\\infty\\left(a_n\\cos\\left(\\frac{2\\pi nx}{P}\\right)+b_n\\sin\\left(\\frac{2\\pi nx}{P}\\right)\\right)").move_to(LOC).scale(0.75).shift(RIGHT),
			TexMobject("a_n=", "\\frac{2}{P}", "\\int_{P}", "s(x)", "\\cos\\left(", "\\frac{2\\pi xn}{P}", "\\right) dx").next_to(LOC, DOWN, buff = 0.5).scale(0.75),
			TexMobject("b_n=\\frac{2}{P}\\int_{P} s(x)\\sin\\left(\\frac{2\\pi xn}{P}\\right) dx").next_to(LOC, DOWN, buff = 1.5).scale(0.75),
			TexMobject("s_N(x)", "=", "\\frac{a_0}{2}", "+", "\\sum_{n=1}^\\infty a_n\\cos\\left(\\frac{2\\pi n(x-5)}{P}\\right)").move_to(LOC).scale(0.75).shift(RIGHT),
			TexMobject("a_n=", "\\frac{1}{5}", "\\int_{-5}^5", "s(x)", "\\cos\\left(", "\\frac{\\pi xn}{5}", "\\right) dx").next_to(LOC, DOWN, buff = 0.5).scale(0.75),
			TexMobject("a_0 = \\frac{37}{6}").next_to(LOC, DOWN, buff = 1.5).scale(0.75),
			TexMobject("s_N(x)", "=", "\\frac{37}{12}", "+", "\\sum_{n=1}^\\infty a_n\\cos\\left(\\frac{2\\pi n(x-5)}{P}\\right)").move_to(LOC).scale(0.75).shift(RIGHT),
			TexMobject("a_n=", "\\frac{1}{5}", "\\int_{-5}^5", "\\left(\\frac{1}{4}x^2+1\\right)", "\\cos\\left(", "\\frac{\\pi xn}{5}", "\\right) dx").next_to(LOC, DOWN, buff = 0.5).scale(0.75),
			TexMobject("a_n = \\frac{25}{\\pi^2n^2}", "\\cos(\\pi n)").next_to(LOC, DOWN, buff = 1.5).scale(0.75),
			TexMobject("a_n = \\frac{25}{\\pi^2n^2}", "(-1)^n").next_to(LOC, DOWN, buff = 1.5).scale(0.75),
			TexMobject("s_N(x)", "=", "\\frac{37}{12}", "+", "\\sum_{n=1}^\\infty", "a_n", "\\cos\\left(\\frac{2\\pi n(x-5)}{P}\\right)").move_to(LOC).scale(0.75).shift(RIGHT),
			TexMobject("s_N(x)", "=", "\\frac{37}{12}", "+", "\\sum_{n=1}^\\infty", "\\frac{(-1)^n\\cdot 25}{\\pi^2n^2}", "\\cos\\left(\\frac{2\\pi n(x-5)}{P}\\right)").move_to(LOC).scale(0.75).shift(RIGHT),
			TexMobject("s_N(x)", "=", "\\frac{37}{12}", "+", "\\sum_{n=1}^", "\\infty", "\\frac{(-1)^n\\cdot 25}{\\pi^2n^2}", "\\cos\\left(\\frac{2\\pi n(x-5)}{P}\\right)").move_to(LOC).scale(0.75).shift(RIGHT),
			TexMobject("s_N(x)", "=", "\\frac{37}{12}", "+", "\\sum_{n=1}^", "k", "\\frac{(-1)^n\\cdot 25}{\\pi^2n^2}", "\\cos\\left(\\frac{2\\pi n(x-5)}{P}\\right)").move_to(LOC).scale(0.75).shift(RIGHT),
		]
		changes = [
			[[(0, 1, 2, 2, 2), (0, 1, 2, 3, 4)]],
			[[(0, 1, 2, 3, 4), (0, 1, 2, 3, 4)]],
			[[(0, 1, 2, 3, 4, 5, 6), (0, 1, 2, 3, 4, 5, 6)]]
		]
		# self.play(Write(eqn[0]))
		# for pre_ind, post_ind in changes[0]:
		# 	self.play(*[Transform(eqn[0][i], eqn[1][j]) for i,j in zip(pre_ind, post_ind)])
		self.play(Write(eqn[1]))
		self.wait(5)
		self.play(Write(eqn[2]))
		self.play(Write(eqn[3]))
		self.wait(5)
		for pre_ind, post_ind in changes[1]:
			self.play(*[Transform(eqn[1][i], eqn[4][j]) for i,j in zip(pre_ind, post_ind)], FadeOut(eqn[3]))
		for pre_ind, post_ind in changes[2]:
			self.play(*[Transform(eqn[2][i], eqn[5][j]) for i,j in zip(pre_ind, post_ind)])
		self.play(FadeOut(graphlabel), FadeIn(shiftedgraphlabel))
		self.wait(5)
		self.play(Transform(eqn[2], eqn[8]))
		self.play(Write(eqn[6]))
		self.wait()
		self.play(Transform(eqn[1], eqn[7]), FadeOut(eqn[6]))
		self.wait(3)
		self.play(Write(eqn[9]))
		self.wait(2)
		self.play(Transform(eqn[9], eqn[10]))
		self.wait(2)
		self.remove(eqn[1])
		self.add(eqn[11])
		self.wait()
		self.play(Transform(eqn[11], eqn[12]))
		self.remove(eqn[11])
		self.add(eqn[13])
		self.play(FadeOut(eqn[9]), FadeOut(eqn[2]), Transform(shiftedgraphlabel, shiftedshiftedgraphlabel))
		self.wait()
		counter = [
			TexMobject("k = ", "?").move_to([5, 0, 0]),
			TexMobject("k = ", "0").move_to([5, 0, 0]),
			TexMobject("k = ", "1").move_to([5, 0, 0]),
			TexMobject("k = ", "2").move_to([5, 0, 0]),
			TexMobject("k = ", "3").move_to([5, 0, 0]),
			TexMobject("k = ", "4").move_to([5, 0, 0]),
			TexMobject("k = ", "5").move_to([5, 0, 0]),
			TexMobject("k = ", "10").move_to([5, 0, 0]),
			TexMobject("k = ", "20").move_to([5, 0, 0]),
			TexMobject("k = ", "50").move_to([5, 0, 0]),
			TexMobject("k = ", "100").move_to([5, 0, 0]),
			TexMobject("k = ", "1000").move_to([5, 0, 0]),
		]
		self.play(Transform(eqn[13], eqn[14]), Write(counter[0]))
		GRAPHCOLOR = "#5406F9"
		grapharr = [
			self.get_graph(fourier0,x_min=0,x_max=10, color = GRAPHCOLOR),
			self.get_graph(fourier1,x_min=0,x_max=10, color = GRAPHCOLOR),
			self.get_graph(fourier2,x_min=0,x_max=10, color = GRAPHCOLOR),
			self.get_graph(fourier3,x_min=0,x_max=10, color = GRAPHCOLOR),
			self.get_graph(fourier4,x_min=0,x_max=10, color = GRAPHCOLOR),
			self.get_graph(fourier5,x_min=0,x_max=10, color = GRAPHCOLOR),
			self.get_graph(fourier10,x_min=0,x_max=10, color = GRAPHCOLOR),
			self.get_graph(fourier20,x_min=0,x_max=10, color = GRAPHCOLOR),
			self.get_graph(fourier50,x_min=0,x_max=10, color = GRAPHCOLOR),
			self.get_graph(fourier100,x_min=0,x_max=10, color = GRAPHCOLOR),
			self.get_graph(fourier1000,x_min=0,x_max=10, color = GRAPHCOLOR),
		]
		for i in range(11):
			self.play(ShowCreation(grapharr[i]), Transform(counter[0], counter[i+1]), run_time = 2)
			self.wait()
			if i != 10:
				self.remove(grapharr[i])
		self.wait()








