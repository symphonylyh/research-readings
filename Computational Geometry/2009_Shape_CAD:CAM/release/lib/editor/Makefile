# Copyright (C) 1995-2008 Massachusetts Institute of Technology, Cambridge, MA
# All rights reserved 

INCL = ../../include 
CFLAGS = $(DEBUG) -I$(INCL)
FFLAGS = $(DEBUG)
EDITOR = ../libeditor.a

OBJS = alpha_shape.o approxcu.o bezierop.o blend_contpts.o \
	blend_cross_link.o blend_curv.o blend_sample.o blend_util.o \
	bound.o bspl_to_bezier.o calc_disc.o check.o compare.o \
	construct.o contour.o convex_hull.o convex_jbspl.o copy.o \
	curv_approx.o decasteljau_curve.o delaunay.o \
        delaunay_to_edge.o deriv.o \
	deviation.o edge.o efg_array.o evaluate.o fair_per_knot.o \
	farray.o find_deviation.o free_earray.o gencyl.o \
	gencyl_util.o get_error.o hull.o intcyl.o \
	integral_offset.o interp_curv.o intersect_curve.o \
	jbspl_util.o knot_rem.o knotv.o knotv_per.o leading_edge.o \
	locate.o loft_integral.o loft_rational.o log.o \
	max_curvature.o nagparam.o normalize.o oslo_per.o partan.o \
	path.o perturb.o prune_facets.o ray.o read.o reflect.o \
	rev_param.o rgb.o ruled.o save_file.o split_surf.o \
	stitch_facet.o subdivide.o surf_offset.o surf_rev.o tab_cyl.o \
	toler_corner.o toler_gencyl.o toler_merge.o toler_region.o \
	transf.o transform.o trim.o trim_contour.o tripack-c.o \
	tripack-f.o vector.o vector_util.o write.o

$(EDITOR): Makefile $(OBJS)
	rm -f $(EDITOR)
	ar ruv $@ $(OBJS)

# force recompile and relink
remake:
	rm -f *.o

# delete backups (*~) and objects (*.o)
clean:
	rm -f *.o *~ \#*





