
"""
	Works well, but not for too many vertices!
	
	Problem:
	If shape is super complex (has many vetrices),
	than cell definition will be too long for MCNP.
	E.g. 1000 vertices of circle is too many.

	Solution:
	* Add option to split whole region into sections.
	This would be done along symetry surfaces.
	Neighbouring cells would share these surfaces.
	This is no problem because of the Cycle class.
	* The user would decide if they want single-cell
	or multicell.
	* Furthermore, the script could export an object
	which is the colection of these cells.
	
	If this is the case, the Offset method should be changed:
	- An offset would not make RegionIn AND RegionOut,
	but instead just new surfaces, each one binding just one cell.
	- The function would instead return a coletion of cells,
	i.e. closed regions.

"""

import numpy as np
#! import openmc
# import matplotlib.pyplot as plt
import cmath
import re


surface_names = [] # list of int
surface_definitions = []



###########################################
#                                         #
#   MCNP SPECIFIC CLASSES AND FUNCTIONS   #
#                                         #
###########################################


class Surfaces:

	def __init__(self):
		global surface_names
		global surface_definitions
		self.surface_definitions  = surface_definitions
		self.surface_names = surface_names

	def FormatSurfaces(self):
		global surface_names
		global surface_definitions
		surface_names = []
		surface_definitions = []

	def __str__(self):
		return '\n'.join(surface_definitions)



class Surface:
	
	def __init__(self, type:str, parameters:list):
		global surface_names
		global surface_definitions
		if len(surface_names) == 0:
			name = 1
		else:
			name = surface_names[-1] + 1
		self.name = name
		surface_names.append(name)
		self.type = type
		self.parameters = parameters
		# return self.name

		surface_definitions.append(str(self))

	def __str__(self):
		parameters = " ".join([str(f) for f in self.parameters])
		return f'{self.name} {self.type} {parameters}'



def Union(regions):
	regions = [str(r) for r in regions]
	definition = '(' + ":".join(regions) + ')'
	return definition

def Intersection(regions):
	regions = [str(r) for r in regions]
	definition = '(' + " ".join(regions) + ')'
	return definition

def is_splitpoint_char(char):
	""" Function for function split_into_lines(),
	which is used by funciton MCNP_format_region(). """
	if char in [' ', '(', ')', ':']:
		return True
	else:
		return False

def split_into_lines(string, max_len):
	""" Function for function MCNP_format_region(),
	which is used by """
	s = string
	list_of_lines = []
	pattern=r'[\+|\-]\d+'

	# check if entity lengths are too long
	margin = 5
	all_entities = re.findall(pattern, string)
	lengths = [len(entity) for entity in all_entities]
	if max(lengths) >= max_len - margin:
		print("Error: Length of entities within margin of max_len.")
		return None

	while len(s) > max_len:
		for i,c in enumerate(s):
			if i+1 > max_len:
				list_of_lines.append(s[:last_possible_splitpoint])
				s = s[last_possible_splitpoint:]
				break
			elif is_splitpoint_char(c):
				last_possible_splitpoint = i
				
	if len(s) != 0:
		# write also the left-overs
		list_of_lines.append(s)
	
	return list_of_lines




###########################################
#                                         #
#            POLYGONTORUS CORE            #
#                                         #
###########################################




class Cycle(object):

	def __init__(self, elements:list):
		self.elements = list(elements)
		self.length = len(elements)

	def __getitem__(self, index:int):
		inside_index = index % self.length
		return self.elements[inside_index]
	
	def RemoveElement(self, index:int):
		if index >= self.length:
			print(f'Warning: Index {index} is outside cycle base and will not be removed.')
		head = self.elements[:index]
		tail = self.elements[index+1:]
		self.elements = head + tail
		self.length = len(self.elements)

	def __str__(self):
		text = f'Cycle object of length {self.length} with elements:'
		for e in self.elements:
			text += f'\n\t{e}'
		return text



class Point:

	def __init__(self, xy, _id):
		self.xy = xy
		self.angle = None
		self.id = _id

	def SetAngle(self, angle):
		self.angle = angle

	def __str__(self):
		return f'Point object:\tid={self.id}\txy={self.xy}\tangle={self.angle} deg'



class Edge:

	def __init__(self, xy1, xy2):
		self.start_point = xy1
		self.end_point = xy2
		self.type = self.GetEdgeType()
		# self.surfaces = None
		self.region_in = None
		self.region_out = None
		self.CalcRegion()

	def GetEdgeType(self):
		x1, y1 = self.start_point
		x2, y2 = self.end_point
		# VERTICALITY
		if x1 < x2: verticality = 'T'
		elif x1 > x2: verticality = 'B'
		else: verticality = ''
		# RADIALITY
		if y1 < y2: radiality = 'I'
		elif y1 > y2: radiality = 'O'
		else: radiality = ''
		return verticality + radiality
	
	###########################################################################
	
	def CalcRegion(self):
		x1, y1 = self.start_point
		x2, y2 = self.end_point
		_t = self.type
		# SLOPED edge
		if _t in ('TI', 'TO', 'BI', 'BO'):
			if 'B' in _t:
				# FLIP IF BOTTOM
				_x1, _y1 = self.end_point
				_x2, _y2 = self.start_point
			else:
				_x1, _y1 = self.start_point
				_x2, _y2 = self.end_point
			k = (_y2-_y1) / (_x2-_x1)
			y0 = _y1 - k * _x1
			t2 = (1/k)**2
			#! cone_surface = openmc.ZCone(x0=0,y0=0, z0=y0, r2=t2)
			cone_surface = Surface('KZ', [y0, t2])
			#! plane_surface = openmc.ZPlane(z0=y0)
			plane_surface = Surface('PZ', [y0])
			# plane region
			if 'T' in _t:
				#! in_plane_region =  -plane_surface
				in_plane_region = f'-{plane_surface.name}'
				#! out_plane_region = +plane_surface
				out_plane_region = f'+{plane_surface.name}'
			elif 'B' in _t:
				#! in_plane_region =  +plane_surface
				in_plane_region = f'+{plane_surface.name}'
				#! out_plane_region = -plane_surface
				out_plane_region = f'-{plane_surface.name}'
			# cone region and total inside region
			if 'I' in _t:
				# SET REGION
				#! self.region_in = openmc.Union([         +cone_surface, in_plane_region ])
				self.region_in = Union([ f'+{cone_surface.name}' , in_plane_region])
				#! self.region_out = openmc.Intersection([ -cone_surface, out_plane_region ])
				self.region_out = Intersection([ f'-{cone_surface.name}', out_plane_region ])
				pass
			elif 'O' in _t:
				# SET REGION
				#! self.region_in = openmc.Intersection([ -cone_surface, in_plane_region ])
				self.region_in = Intersection([ f'-{cone_surface.name}', in_plane_region ])
				#! self.region_out = openmc.Union([       +cone_surface, out_plane_region ])
				self.region_out = Union([       f'+{cone_surface.name}', out_plane_region ])
				pass
		# HORIZONTAL edge
		elif _t in ('T', 'B'):
			y0 = y1
			#! plane_surface = openmc.ZPlane(z0=y0) #  f'plane {y0}'
			plane_surface = Surface('PZ', [y0]) #  f'plane {y0}'
			if _t == 'B':
				#! self.region_in =  +plane_surface
				self.region_in =  f'+{plane_surface.name}'
				#! self.region_out = -plane_surface
				self.region_out = f'-{plane_surface.name}'
			elif _t == 'T':
				#! self.region_in =  -plane_surface
				self.region_in =  f'-{plane_surface.name}'
				#! self.region_out = +plane_surface
				self.region_out = f'+{plane_surface.name}'
		 # VERTICAL edge
		elif _t in ('I', 'O'):
			x0 = x1
			#! cylinder_surface = openmc.ZCylinder(x0=0,y0=0, r=x0)
			cylinder_surface = Surface('CZ', [x0])
			if _t == 'I':
				#! self.region_in =  +cylinder_surface
				self.region_in =  f'+{cylinder_surface.name}'
				#! self.region_out = -cylinder_surface
				self.region_out = f'-{cylinder_surface.name}'
			elif _t == 'O':
				#! self.region_in =  -cylinder_surface
				self.region_in =  f'-{cylinder_surface.name}'
				#! self.region_out = +cylinder_surface
				self.region_out = f'+{cylinder_surface.name}'
		else:
			print(f'Error: "{self.type}" is unknown edge type!')
	
	###########################################################################
	
	def __str__(self):
		text = f'Edge object:'
		text += f'\n\tFrom {self.start_point} to {self.end_point}'
		text += f'\n\tType = {self.type}'
		text += f'\n\tRegion inside = {self.region_in}'
		text += f'\n\tRegion inside = {self.region_out}'
		return text



class Polygon:

	def __init__(self, xy_points, remove_nonconvex=True):
		# self.points = xy_points
		points = [Point(xy=xy, _id=i) for i,xy in enumerate(xy_points)]
		self.points = Cycle(points)
		self.CalcPointAngles()
		self.original_xy_points = xy_points
		if remove_nonconvex:
			self.RemoveNonConvexPoints()
		self.edges = self.GetEdges()

	def GetEdges(self):
		Edges = []
		for i in range(self.points.length):
			xy1 = self.points[i].xy
			xy2 = self.points[i+1].xy
			Edges.append(Edge(xy1, xy2))
		return Edges
	
	def GetAngleToNeighbouringPoints(self, point_index:int):
		""" Returns angle to previous point and angle to next point.
		Angle is in the sense of a point in the complex plane with the origin set to point_index point. """
		i = point_index
		point_p = self.points[i-1] # previous point
		z_p = complex(point_p.xy[0], point_p.xy[1])
		point_c = self.points[i] # current point
		z_c = complex(point_c.xy[0], point_c.xy[1])
		point_n = self.points[i+1] # next point
		z_n = complex(point_n.xy[0], point_n.xy[1])
		# set origin to current point
		z_PP = z_p - z_c
		z_NN = z_n - z_c
		angle_p = np.angle(z_PP, deg=True)
		angle_n = np.angle(z_NN, deg=True)
		return (angle_p, angle_n)
	
	def CalcPointAngles(self):
		for i in range(self.points.length):
			angle_p, angle_n = self.GetAngleToNeighbouringPoints(point_index=i)
			if angle_p > angle_n:
				angle = (360 - angle_p) + angle_n
			else:
				angle = angle_n - angle_p
			self.points[i].SetAngle(angle)

	def RemoveNonConvexPoints(self):
		""" Method to remove all co-linear or concave points from shell.
		The concavity of the points is determined by the inside angle of the two neighbouring edges.
		If this angle is 180 degrees or more the point is non-covenx and it should be removed. """
		i = 0
		n_points = self.points.length
		already_removed_some = False
		while i < n_points:
			if self.points.elements[i].angle >= 180:
				if not already_removed_some:
					print('Removing non-convex points:')
				print(f'\t{self.points.elements[i]}')
				self.points.RemoveElement(i)
				i = 0 # go from start again
				n_points = self.points.length
				self.CalcPointAngles() # recalculate point angles since we're changing them !
				already_removed_some = True
			else:
				i += 1 # check next point
			

	def Offset(self, d:float):
		""" Function returns new Polygon object with offset points.
		Points are offset in the direction of the simetral axis.
		The simetral axis is determined by the neighbouring edges. """
		new_points_xy = []
		for i in range(self.points.length):
			point_c = self.points[i] # current point
			angle_c = point_c.angle
			ac_rad = np.deg2rad(angle_c)
			z_c = complex(point_c.xy[0], point_c.xy[1])
			angle_p, angle_n = self.GetAngleToNeighbouringPoints(point_index=i)
			mean_angle = (angle_n + angle_p) / 2
			if angle_p > angle_n:
				offset_angle = mean_angle
			else:
				offset_angle = mean_angle - 180 # deg
			oa_rad = np.deg2rad(offset_angle)
			offset = d / np.sin(ac_rad / 2)
			z_offset = cmath.rect(offset, oa_rad)
			new_z = z_c + z_offset
			new_points_xy.append([new_z.real, new_z.imag])
		return new_points_xy
	
	def PlotPolygon(self, mpl_axes, **plot_kwargs):
		pt_arr = np.array([point.xy for point in self.points.elements])
		pt_arr = np.concatenate((pt_arr, [pt_arr[0]]))
		x = pt_arr[:,0]
		y = pt_arr[:,1]
		mpl_axes.plot(x,y, **plot_kwargs)
	
	def __str__(self):
		text = f'Polygon object with {self.points.length} perimeter points:'
		for p in self.points.elements:
			text += f'\n\txy={p.xy}\tangle={p.angle} deg'
		return text



class PolygonTorus:

	def __init__(self, rz_points):
		self.polygon = Polygon(rz_points)
		self.region_in = self.GetRegionInside()
		""" self.region_in = '\n     '.join(split_into_lines(_region_in,
												 max_len=74)) """
		self.region_out = self.GetRegionOutside()
		""" self.region_out = '\n     '.join(split_into_lines(_region_out,
												 max_len=74)) """
	
	
	def GetRegionInside(self):
		edges_inside_regions = [e.region_in for e in self.polygon.edges]
		# regions_str = '\n\t'.join(edges_inside_regions)
		#! return openmc.Intersection(edges_inside_regions) # f'INTERSECTION(\n\t{regions_str}\n)'# 
		return Intersection(edges_inside_regions) # f'INTERSECTION(\n\t{regions_str}\n)'# 
	
	def GetRegionOutside(self):
		edges_outside_regions = [e.region_out for e in self.polygon.edges]
		# regions_str = '\n\t'.join(edges_outside_regions)
		#! return openmc.Union(edges_outside_regions) # f'UNION(\n\t{regions_str}\n)'# 	
		return Union(edges_outside_regions) # f'UNION(\n\t{regions_str}\n)'# 	
	
	def Offset(self, d:float):
		""" Function returns new PolygonTorus object with offset shell.
		This function calls the Polygon.Offset(d) function. """
		new_xy_points = self.polygon.Offset(d)
		new_polygon_torus = PolygonTorus(new_xy_points)
		current_PT_n_points = self.polygon.points.length
		offset_PT_n_points = new_polygon_torus.polygon.points.length
		if offset_PT_n_points < current_PT_n_points:
			print(f'Warrning: Offset PolygonTorus has {current_PT_n_points-offset_PT_n_points} less points than original PolygonTorus. This is a result of setting a large negative offset.')
		return new_polygon_torus
	
	def __str__(self):
		text = f'PolygonTorus object:'
		text += f'\n\tPoints: {self.polygon.points.elements}'
		text += f'\n\tEdge types are: {[e.type for e in self.polygon.edges]}'
		text += f'\n\tRegion inside is: {self.region_in}'
		text += f'\n\tRegion outside is: {self.region_out}'
		return text

	def GetSurfaces(self):
		""" Returns list of surfaces created in session. """
		return surface_definitions
	