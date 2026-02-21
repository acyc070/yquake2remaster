/*
 * Copyright (C) 1997-2001 Id Software, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 *
 */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include "header/local.h"

#define AFFINE_SPANLET_SIZE      16
#define AFFINE_SPANLET_SIZE_BITS 4

typedef struct
{
	pixel_t		*pbase, *pdest;
	zvalue_t	*pz;
	int		s, t;
	int		sstep, tstep;
	int		izi, izistep, izistep_times_2;
	int		spancount;
	unsigned	u, v;
} spanletvars_t;

static spanletvars_t	s_spanletvars;

static int	r_polyblendcolor;

polydesc_t	r_polydesc;

msurface_t	*r_alpha_surfaces;

static int	clip_current;
vec5_t		r_clip_verts[2][MAXWORKINGVERTS+2];
static emitpoint_t	outverts[MAXWORKINGVERTS+3];

static int	s_minindex, s_maxindex;

// PS1-style affine texture mapping coefficients
static float	s_A, s_B, s_C;
static float	t_A, t_B, t_C;

static void R_DrawPoly(int iswater, espan_t *spans);

/*
** R_DrawSpanletAffineOpaque - PS1 style affine texture mapping
*/
static void
R_DrawSpanletAffineOpaque(const int *r_turb_turb)
{
	do
	{
		unsigned btemp;
		unsigned ts, tt;

		// Affine mapping - just use interpolated s,t directly
		ts = s_spanletvars.s >> SHIFT16XYZ;
		tt = s_spanletvars.t >> SHIFT16XYZ;

		// Clamp to texture bounds for PS1-style wrapping
		ts &= (cachewidth - 1);
		tt &= (cacheheight - 1);

		btemp = *(s_spanletvars.pbase + (ts) + (tt) * cachewidth);
		if (btemp != TRANSPARENT_COLOR)
		{
			if (*s_spanletvars.pz <= (s_spanletvars.izi >> SHIFT16XYZ))
			{
				*s_spanletvars.pz    = s_spanletvars.izi >> SHIFT16XYZ;
				*s_spanletvars.pdest = btemp;
			}
		}

		s_spanletvars.izi += s_spanletvars.izistep;
		s_spanletvars.pdest++;
		s_spanletvars.pz++;
		s_spanletvars.s += s_spanletvars.sstep;
		s_spanletvars.t += s_spanletvars.tstep;
	} while (--s_spanletvars.spancount > 0);
}

/*
** R_DrawSpanletAffine33 - PS1 style with alpha blending
*/
static void
R_DrawSpanletAffine33(const int *r_turb_turb)
{
	do
	{
		unsigned btemp;
		unsigned ts, tt;

		ts = s_spanletvars.s >> SHIFT16XYZ;
		tt = s_spanletvars.t >> SHIFT16XYZ;

		// PS1-style texture clamping/wrapping
		ts &= (cachewidth - 1);
		tt &= (cacheheight - 1);

		btemp = *(s_spanletvars.pbase + (ts) + (tt) * cachewidth);

		if (btemp != TRANSPARENT_COLOR)
		{
			if (*s_spanletvars.pz <= (s_spanletvars.izi >> SHIFT16XYZ))
			{
				*s_spanletvars.pdest = vid_alphamap[btemp + *s_spanletvars.pdest * 256];
			}
		}

		s_spanletvars.izi += s_spanletvars.izistep;
		s_spanletvars.pdest++;
		s_spanletvars.pz++;
		s_spanletvars.s += s_spanletvars.sstep;
		s_spanletvars.t += s_spanletvars.tstep;
	} while (--s_spanletvars.spancount > 0);
}

/*
** R_DrawSpanletAffine66 - PS1 style with alpha blending (66 mode)
*/
static void
R_DrawSpanletAffine66(const int *r_turb_turb)
{
	do
	{
		unsigned btemp;
		unsigned ts, tt;

		ts = s_spanletvars.s >> SHIFT16XYZ;
		tt = s_spanletvars.t >> SHIFT16XYZ;

		// PS1-style texture clamping/wrapping
		ts &= (cachewidth - 1);
		tt &= (cacheheight - 1);

		btemp = *(s_spanletvars.pbase + (ts) + (tt) * cachewidth);

		if (btemp != TRANSPARENT_COLOR)
		{
			if (*s_spanletvars.pz <= (s_spanletvars.izi >> SHIFT16XYZ))
			{
				*s_spanletvars.pdest = vid_alphamap[btemp * 256 + *s_spanletvars.pdest];
			}
		}

		s_spanletvars.izi += s_spanletvars.izistep;
		s_spanletvars.pdest++;
		s_spanletvars.pz++;
		s_spanletvars.s += s_spanletvars.sstep;
		s_spanletvars.t += s_spanletvars.tstep;
	} while (--s_spanletvars.spancount > 0);
}

/*
** R_DrawSpanletAffineConstant33 - PS1 style flat shading
*/
static void
R_DrawSpanletAffineConstant33(const int *r_turb_turb)
{
	do
	{
		if (*s_spanletvars.pz <= (s_spanletvars.izi >> SHIFT16XYZ))
		{
			*s_spanletvars.pdest = vid_alphamap[r_polyblendcolor + *s_spanletvars.pdest * 256];
		}

		s_spanletvars.izi += s_spanletvars.izistep;
		s_spanletvars.pdest++;
		s_spanletvars.pz++;
	} while (--s_spanletvars.spancount > 0);
}

/*
** R_DrawSpanletAffine33Stipple - PS1 style stippled transparency
*/
static void
R_DrawSpanletAffine33Stipple(const int *r_turb_turb)
{
	pixel_t		*pdest	= s_spanletvars.pdest;
	zvalue_t	*pz	= s_spanletvars.pz;
	zvalue_t	izi	= s_spanletvars.izi;

	if ( r_polydesc.stipple_parity ^ ( s_spanletvars.v & 1 ) )
	{
		s_spanletvars.pdest += s_spanletvars.spancount;
		s_spanletvars.pz    += s_spanletvars.spancount;

		if ( s_spanletvars.spancount == AFFINE_SPANLET_SIZE )
			s_spanletvars.izi += s_spanletvars.izistep << AFFINE_SPANLET_SIZE_BITS;
		else
			s_spanletvars.izi += s_spanletvars.izistep * s_spanletvars.izistep;

		if ( r_polydesc.stipple_parity ^ ( s_spanletvars.u & 1 ) )
		{
			izi += s_spanletvars.izistep;
			s_spanletvars.s   += s_spanletvars.sstep;
			s_spanletvars.t   += s_spanletvars.tstep;

			pdest++;
			pz++;
			s_spanletvars.spancount--;
		}

		s_spanletvars.sstep *= 2;
		s_spanletvars.tstep *= 2;

		while ( s_spanletvars.spancount > 0 )
		{
			unsigned btemp;
			unsigned s = s_spanletvars.s >> SHIFT16XYZ;
			unsigned t = s_spanletvars.t >> SHIFT16XYZ;

			// PS1-style texture clamping/wrapping
			s &= (cachewidth - 1);
			t &= (cacheheight - 1);

			btemp = *( s_spanletvars.pbase + ( s ) + ( t * cachewidth ) );

			if ( btemp != TRANSPARENT_COLOR )
			{
				if ( *pz <= ( izi >> SHIFT16XYZ ) )
					*pdest = btemp;
			}

			izi               += s_spanletvars.izistep_times_2;
			s_spanletvars.s   += s_spanletvars.sstep;
			s_spanletvars.t   += s_spanletvars.tstep;

			pdest += 2;
			pz    += 2;

			s_spanletvars.spancount -= 2;
		}
	}
}

/*
** R_DrawSpanletAffine66Stipple - PS1 style stippled transparency (66 mode)
*/
static void
R_DrawSpanletAffine66Stipple(const int *r_turb_turb)
{
	unsigned	btemp;
	pixel_t		*pdest = s_spanletvars.pdest;
	zvalue_t	*pz = s_spanletvars.pz;
	zvalue_t	izi = s_spanletvars.izi;

	s_spanletvars.pdest += s_spanletvars.spancount;
	s_spanletvars.pz    += s_spanletvars.spancount;

	if ( s_spanletvars.spancount == AFFINE_SPANLET_SIZE )
		s_spanletvars.izi += s_spanletvars.izistep << AFFINE_SPANLET_SIZE_BITS;
	else
		s_spanletvars.izi += s_spanletvars.izistep * s_spanletvars.izistep;

	if ( r_polydesc.stipple_parity ^ ( s_spanletvars.v & 1 ) )
	{
		if ( r_polydesc.stipple_parity ^ ( s_spanletvars.u & 1 ) )
		{
			izi += s_spanletvars.izistep;
			s_spanletvars.s += s_spanletvars.sstep;
			s_spanletvars.t += s_spanletvars.tstep;

			pdest++;
			pz++;
			s_spanletvars.spancount--;
		}

		s_spanletvars.sstep *= 2;
		s_spanletvars.tstep *= 2;

		while ( s_spanletvars.spancount > 0 )
		{
			unsigned s = s_spanletvars.s >> SHIFT16XYZ;
			unsigned t = s_spanletvars.t >> SHIFT16XYZ;

			// PS1-style texture clamping/wrapping
			s &= (cachewidth - 1);
			t &= (cacheheight - 1);

			btemp = *( s_spanletvars.pbase + ( s ) + ( t * cachewidth ) );

			if ( btemp != TRANSPARENT_COLOR )
			{
				if ( *pz <= ( izi >> SHIFT16XYZ ) )
					*pdest = btemp;
			}

			izi             += s_spanletvars.izistep_times_2;
			s_spanletvars.s += s_spanletvars.sstep;
			s_spanletvars.t += s_spanletvars.tstep;

			pdest += 2;
			pz    += 2;

			s_spanletvars.spancount -= 2;
		}
	}
	else
	{
		while ( s_spanletvars.spancount > 0 )
		{
			unsigned s = s_spanletvars.s >> SHIFT16XYZ;
			unsigned t = s_spanletvars.t >> SHIFT16XYZ;

			// PS1-style texture clamping/wrapping
			s &= (cachewidth - 1);
			t &= (cacheheight - 1);

			btemp = *( s_spanletvars.pbase + ( s ) + ( t * cachewidth ) );

			if ( btemp != TRANSPARENT_COLOR )
			{
				if ( *pz <= ( izi >> SHIFT16XYZ ) )
					*pdest = btemp;
			}

			izi             += s_spanletvars.izistep;
			s_spanletvars.s += s_spanletvars.sstep;
			s_spanletvars.t += s_spanletvars.tstep;

			pdest++;
			pz++;

			s_spanletvars.spancount--;
		}
	}
}

/*
** R_ClipPolyFace
**
** Clips the winding at clip_verts[clip_current] and changes clip_current
** Throws out the back side
*/
static int
R_ClipPolyFace(int nump, const clipplane_t *pclipplane)
{
	float dists[MAXWORKINGVERTS+3] = {0};
	const float *vert2, *pclipnormal;
	float *in, *instep, *outstep;
	float frac, clipdist;
	int i, outcount;

	clipdist = pclipplane->dist;
	pclipnormal = pclipplane->normal;

	// calc dists
	if (clip_current)
	{
		in = r_clip_verts[1][0];
		outstep = r_clip_verts[0][0];
		clip_current = 0;
	}
	else
	{
		in = r_clip_verts[0][0];
		outstep = r_clip_verts[1][0];
		clip_current = 1;
	}

	instep = in;
	for (i=0 ; i<nump ; i++, instep += sizeof (vec5_t) / sizeof (vec_t))
	{
		dists[i] = DotProduct (instep, pclipnormal) - clipdist;
	}

	/* handle wraparound case */
	dists[nump] = dists[0];
	memmove(instep, in, sizeof (vec5_t));

	/* clip the winding */
	instep = in;
	outcount = 0;

	for (i=0 ; i<nump ; i++, instep += sizeof (vec5_t) / sizeof (vec_t))
	{
		if (dists[i] >= 0)
		{
			memcpy (outstep, instep, sizeof (vec5_t));
			outstep += sizeof (vec5_t) / sizeof (vec_t);
			outcount++;
		}

		if (dists[i] == 0 || dists[i+1] == 0)
			continue;

		if ( (dists[i] > 0) == (dists[i+1] > 0) )
			continue;

		// split it into a new vertex
		frac = dists[i] / (dists[i] - dists[i+1]);

		vert2 = instep + sizeof (vec5_t) / sizeof (vec_t);

		outstep[0] = instep[0] + frac*(vert2[0] - instep[0]);
		outstep[1] = instep[1] + frac*(vert2[1] - instep[1]);
		outstep[2] = instep[2] + frac*(vert2[2] - instep[2]);
		outstep[3] = instep[3] + frac*(vert2[3] - instep[3]);
		outstep[4] = instep[4] + frac*(vert2[4] - instep[4]);

		outstep += sizeof (vec5_t) / sizeof (vec_t);
		outcount++;
	}

	return outcount;
}

/*
** R_PolygonDrawSpans - PS1 style affine texture mapping
*/
static void
R_PolygonDrawSpans(espan_t *pspan, int iswater, float d_ziorigin, float d_zistepu, float d_zistepv)
{
	float	du, dv;
	int	*r_turb_turb;

	s_spanletvars.pbase = cacheblock;

	if ( iswater & SURF_WARP)
		r_turb_turb = sintable + ((int)(r_newrefdef.time*SPEED)&(CYCLE-1));
	else
		r_turb_turb = blanktable;

	// PS1 style - use fixed steps for affine mapping
	s_spanletvars.izistep = (int)(d_zistepu * 0x8000 * SHIFT16XYZ_MULT);
	s_spanletvars.izistep_times_2 = s_spanletvars.izistep * 2;

	s_spanletvars.pz = 0;

	do
	{
		int	count;

		s_spanletvars.pdest   = d_viewbuffer + (vid_buffer_width * pspan->v) + pspan->u;
		s_spanletvars.pz      = d_pzbuffer + (vid_buffer_width * pspan->v) + pspan->u;
		s_spanletvars.u       = pspan->u;
		s_spanletvars.v       = pspan->v;
		count = pspan->count;

		if (count > 0)
		{
			// transparent spans damage z buffer
			VID_DamageZBuffer(pspan->u, pspan->v);
			VID_DamageZBuffer(pspan->u + count, pspan->v);

			// PS1 style - calculate initial s,t using affine equations
			du = (float)pspan->u;
			dv = (float)pspan->v;

			// Affine texture coordinate calculation (linear in screen space)
			s_spanletvars.s = (int)((s_A * du + s_B * dv + s_C) * SHIFT16XYZ_MULT);
			s_spanletvars.t = (int)((t_A * du + t_B * dv + t_C) * SHIFT16XYZ_MULT);

			// Keep perspective-correct depth for Z-buffer
			float zi = d_ziorigin + dv * d_zistepv + du * d_zistepu;
			s_spanletvars.izi = (int)(zi * 0x8000 * SHIFT16XYZ_MULT);

			// PS1 style - constant steps per pixel (affine)
			int s_step_pixel = (int)(s_A * SHIFT16XYZ_MULT);
			int t_step_pixel = (int)(t_B * SHIFT16XYZ_MULT);
			
			// For horizontal spans, s changes more than t
			// Adjust for better PS1-style warping
			s_spanletvars.sstep = s_step_pixel;
			s_spanletvars.tstep = t_step_pixel;

			do
			{
				if (count >= AFFINE_SPANLET_SIZE)
					s_spanletvars.spancount = AFFINE_SPANLET_SIZE;
				else
					s_spanletvars.spancount = count;

				count -= s_spanletvars.spancount;

				// PS1 style - just use the same steps for the whole spanlet
				// No perspective correction at spanlet boundaries

				// Call the appropriate spanlet drawing function
				r_polydesc.drawspanlet(r_turb_turb);

				// Update for next spanlet using affine steps
				s_spanletvars.s += s_step_pixel * s_spanletvars.spancount;
				s_spanletvars.t += t_step_pixel * s_spanletvars.spancount;

			} while (count > 0);
		}

		pspan++;

	} while (pspan->count != INT_MIN);
}

/*
** R_PolygonCalculateAffineCoefficients
**
** Calculates affine mapping coefficients for PS1-style texture mapping
** Solves for s = A*u + B*v + C using three vertices
*/
static void
R_PolygonCalculateAffineCoefficients(const emitpoint_t *pverts)
{
	// Use first three vertices to solve for affine coefficients
	float u0 = pverts[0].u;
	float v0 = pverts[0].v;
	float s0 = pverts[0].s;
	float t0 = pverts[0].t;
	
	float u1 = pverts[1].u;
	float v1 = pverts[1].v;
	float s1 = pverts[1].s;
	float t1 = pverts[1].t;
	
	float u2 = pverts[2].u;
	float v2 = pverts[2].v;
	float s2 = pverts[2].s;
	float t2 = pverts[2].t;

	// Calculate denominator for plane equation
	float denom = (u1 - u0) * (v2 - v0) - (u2 - u0) * (v1 - v0);
	
	if (fabs(denom) < 0.0001f)
	{
		// Degenerate triangle - use fallback
		s_A = 0; s_B = 0; s_C = s0;
		t_A = 0; t_B = 0; t_C = t0;
		return;
	}

	// Solve for s coefficients
	s_A = ((s1 - s0) * (v2 - v0) - (s2 - s0) * (v1 - v0)) / denom;
	s_B = ((u1 - u0) * (s2 - s0) - (u2 - u0) * (s1 - s0)) / denom;
	s_C = s0 - s_A * u0 - s_B * v0;

	// Solve for t coefficients
	t_A = ((t1 - t0) * (v2 - v0) - (t2 - t0) * (v1 - v0)) / denom;
	t_B = ((u1 - u0) * (t2 - t0) - (u2 - u0) * (t1 - t0)) / denom;
	t_C = t0 - t_A * u0 - t_B * v0;
}

/*
** R_PolygonScanLeftEdge
**
** Goes through the polygon and scans the left edge, filling in
** screen coordinate data for the spans
*/
static void
R_PolygonScanLeftEdge (espan_t *s_polygon_spans)
{
	const emitpoint_t *pvert, *pnext;
	float du, dv, vtop, u_step;
	int i, lmaxindex;
	espan_t *pspan;

	pspan = s_polygon_spans;
	i = s_minindex;
	if (i == 0)
		i = r_polydesc.nump;

	lmaxindex = s_maxindex;
	if (lmaxindex == 0)
		lmaxindex = r_polydesc.nump;

	vtop = ceil (r_polydesc.pverts[i].v);

	do
	{
		float vbottom;

		pvert = &r_polydesc.pverts[i];
		pnext = pvert - 1;

		vbottom = ceil (pnext->v);

		if (vtop < vbottom)
		{
			int v, u, istep, itop, ibottom;

			du = pnext->u - pvert->u;
			dv = pnext->v - pvert->v;

			u_step = (du * SHIFT16XYZ_MULT) / dv;
			// adjust u to ceil the integer portion
			u = (int)((pvert->u * SHIFT16XYZ_MULT) + u_step * (vtop - pvert->v)) +
				(SHIFT16XYZ_MULT - 1);
			itop = (int)vtop;
			ibottom = (int)vbottom;
			istep = (int)u_step;

			for (v=itop ; v<ibottom ; v++)
			{
				pspan->u = u >> SHIFT16XYZ;
				pspan->v = v;
				u += istep;
				pspan++;
			}
		}

		vtop = vbottom;

		i--;
		if (i == 0)
			i = r_polydesc.nump;

	} while (i != lmaxindex);

	pspan->count = INT_MIN;	// mark the end of the span list
}

/*
** R_PolygonScanRightEdge
**
** Goes through the polygon and scans the right edge, filling in
** count values.
*/
static void
R_PolygonScanRightEdge(espan_t *s_polygon_spans)
{
	float du, dv, vtop, u_step, uvert, unext, vvert;
	const emitpoint_t *pnext;
	emitpoint_t *pvert;
	espan_t *pspan;
	int i;

	pspan = s_polygon_spans;
	i = s_minindex;

	vvert = r_polydesc.pverts[i].v;
	if (vvert < r_refdef.fvrecty_adj)
		vvert = r_refdef.fvrecty_adj;
	if (vvert > r_refdef.fvrectbottom_adj)
		vvert = r_refdef.fvrectbottom_adj;

	vtop = ceil (vvert);

	do
	{
		float vbottom, vnext;

		pvert = &r_polydesc.pverts[i];
		pnext = pvert + 1;

		vnext = pnext->v;
		if (vnext < r_refdef.fvrecty_adj)
			vnext = r_refdef.fvrecty_adj;
		if (vnext > r_refdef.fvrectbottom_adj)
			vnext = r_refdef.fvrectbottom_adj;

		vbottom = ceil (vnext);

		if (vtop < vbottom)
		{
			int v, u, istep, itop, ibottom;

			uvert = pvert->u;
			if (uvert < r_refdef.fvrectx_adj)
				uvert = r_refdef.fvrectx_adj;
			if (uvert > r_refdef.fvrectright_adj)
				uvert = r_refdef.fvrectright_adj;

			unext = pnext->u;
			if (unext < r_refdef.fvrectx_adj)
				unext = r_refdef.fvrectx_adj;
			if (unext > r_refdef.fvrectright_adj)
				unext = r_refdef.fvrectright_adj;

			du = unext - uvert;
			dv = vnext - vvert;
			u_step = (du * SHIFT16XYZ_MULT) / dv;
			// adjust u to ceil the integer portion
			u = (int)((uvert * SHIFT16XYZ_MULT) + u_step * (vtop - vvert)) +
					(SHIFT16XYZ_MULT - 1);
			itop = (int)vtop;
			ibottom = (int)vbottom;
			istep = (int)u_step;

			for (v=itop ; v<ibottom ; v++)
			{
				pspan->count = (u >> SHIFT16XYZ) - pspan->u;
				u += istep;
				pspan++;
			}
		}

		vtop = vbottom;
		vvert = vnext;

		i++;
		if (i == r_polydesc.nump)
			i = 0;

	} while (i != s_maxindex);

	pspan->count = INT_MIN;	// mark the end of the span list
}

/*
** R_ClipAndDrawPoly - PS1 style affine mapping
*/
void
R_ClipAndDrawPoly(float alpha, int isturbulent, qboolean textured)
{
	vec_t		*pv;
	int		i, nump;
	vec3_t		transformed, local;

	if (!textured)
	{
		r_polydesc.drawspanlet = R_DrawSpanletAffineConstant33;
	}
	else
	{
		/*
		** choose the correct spanlet routine based on alpha
		** All using PS1-style affine mapping
		*/
		if ( alpha == 1 )
		{
			r_polydesc.drawspanlet = R_DrawSpanletAffineOpaque;
		}
		else
		{
			if ( sw_stipplealpha->value )
			{
				if ( alpha > 0.33 )
					r_polydesc.drawspanlet = R_DrawSpanletAffine66Stipple;
				else
					r_polydesc.drawspanlet = R_DrawSpanletAffine33Stipple;
			}
			else
			{
				if ( alpha > 0.33 )
					r_polydesc.drawspanlet = R_DrawSpanletAffine66;
				else
					r_polydesc.drawspanlet = R_DrawSpanletAffine33;
			}
		}
	}

	// clip to the frustum in worldspace
	nump = r_polydesc.nump;
	clip_current = 0;

	for (i=0 ; i<4 ; i++)
	{
		nump = R_ClipPolyFace (nump, &view_clipplanes[i]);
		if (nump < 3)
			return;
		if (nump > MAXWORKINGVERTS)
		{
			Com_Error(ERR_DROP, "%s: too many points: %d", __func__, nump);
			return;
		}
	}

	// transform vertices into viewspace and project
	pv = &r_clip_verts[clip_current][0][0];

	for (i=0 ; i<nump ; i++)
	{
		float scale;
		emitpoint_t *pout;

		VectorSubtract (pv, r_origin, local);
		TransformVector (local, transformed);

		if (transformed[2] < NEAR_CLIP)
			transformed[2] = NEAR_CLIP;

		pout = &outverts[i];
		pout->zi = 1.0 / transformed[2];

		pout->s = pv[3];
		pout->t = pv[4];

		scale = xscale * pout->zi;
		pout->u = (xcenter + scale * transformed[0]);

		scale = yscale * pout->zi;
		pout->v = (ycenter - scale * transformed[1]);

		pv += sizeof (vec5_t) / sizeof (vec_t);
	}

	// draw it
	r_polydesc.nump = nump;
	r_polydesc.pverts = outverts;

	R_DrawPoly(isturbulent, vid_polygon_spans);
}

/*
** R_BuildPolygonFromSurface
*/
static void
R_BuildPolygonFromSurface(const entity_t *currententity, const model_t *currentmodel, msurface_t *fa)
{
	medge_t *pedges, *r_pedge;
	float tmins[2] = { 0, 0 };
	int i, lnumverts;
	const float *vec;
	vec5_t *pverts;

	r_polydesc.nump = 0;

	/* reconstruct the polygon */
	pedges = currentmodel->edges;
	lnumverts = fa->numedges;

	pverts = r_clip_verts[0];

	for (i=0 ; i<lnumverts ; i++)
	{
		int lindex;

		lindex = currentmodel->surfedges[fa->firstedge + i];

		if (lindex > 0)
		{
			r_pedge = &pedges[lindex];
			vec = currentmodel->vertexes[r_pedge->v[0]].position;
		}
		else
		{
			r_pedge = &pedges[-lindex];
			vec = currentmodel->vertexes[r_pedge->v[1]].position;
		}

		VectorCopy (vec, pverts[i] );
	}

	VectorCopy( fa->texinfo->vecs[0], r_polydesc.vright );
	VectorCopy( fa->texinfo->vecs[1], r_polydesc.vup );
	VectorCopy( fa->plane->normal, r_polydesc.vpn );
	VectorCopy( r_origin, r_polydesc.viewer_position );

	if ( fa->flags & SURF_PLANEBACK )
	{
		VectorSubtract( vec3_origin, r_polydesc.vpn, r_polydesc.vpn );
	}

	if ( fa->texinfo->flags & (SURF_WARP | SURF_SCROLL) )
	{
		r_polydesc.pixels       = fa->texinfo->image->pixels[0];
		r_polydesc.pixel_width  = fa->texinfo->image->width;
		r_polydesc.pixel_height = fa->texinfo->image->height;
	}
	else
	{
		surfcache_t *scache;

		scache = D_CacheSurface(currententity, fa, 0);

		r_polydesc.pixels       = scache->data;
		r_polydesc.pixel_width  = scache->width;
		r_polydesc.pixel_height = scache->height;

		tmins[0] = fa->texturemins[0];
		tmins[1] = fa->texturemins[1];
	}

	r_polydesc.dist = DotProduct( r_polydesc.vpn, pverts[0] );

	r_polydesc.s_offset = fa->texinfo->vecs[0][3] - tmins[0];
	r_polydesc.t_offset = fa->texinfo->vecs[1][3] - tmins[1];

	// scrolling texture addition
	if (fa->texinfo->flags & SURF_SCROLL)
	{
		float sscroll, tscroll;

		R_FlowingScroll(&r_newrefdef, fa->texinfo->flags, &sscroll, &tscroll);

		r_polydesc.s_offset += 2 * sscroll;
		r_polydesc.t_offset += 2 * tscroll;
	}

	r_polydesc.nump = lnumverts;
}

/*
** R_PolygonCalculateGradients - Keep for Z-buffer only
*/
static void
R_PolygonCalculateGradients (float *p_ziorigin, float *p_zistepu, float *p_zistepv)
{
	vec3_t	p_normal;
	float	distinv;
	float	d_ziorigin, d_zistepu, d_zistepv;

	TransformVector (r_polydesc.vpn, p_normal);

	distinv = 1.0 / (-(DotProduct (r_polydesc.viewer_position, r_polydesc.vpn)) + r_polydesc.dist );

	// Only calculate Z gradients - we don't need s/t gradients anymore
	d_zistepu =   p_normal[0] * xscaleinv * distinv;
	d_zistepv =  -p_normal[1] * yscaleinv * distinv;
	d_ziorigin =  p_normal[2] * distinv - xcenter * d_zistepu - ycenter * d_zistepv;

	*p_zistepu = d_zistepu;
	*p_zistepv = d_zistepv;
	*p_ziorigin = d_ziorigin;
	
	// Texture dimensions for clamping
	bbextents = r_polydesc.pixel_width;
	bbextentt = r_polydesc.pixel_height;
}

/*
** R_DrawPoly - PS1 style affine mapping
**
** Polygon drawing function. Uses the polygon described in r_polydesc
** to calculate edges and gradients, then renders the resultant spans.
*/
static void
R_DrawPoly(int iswater, espan_t *spans)
{
	int		i, nump;
	float		ymin, ymax;
	emitpoint_t	*pverts;
	float	d_ziorigin, d_zistepu, d_zistepv;

	// find the top and bottom vertices, and make sure there's at least one scan to
	// draw
	ymin = (float)INT_MAX;
	ymax = (float)INT_MIN;
	pverts = r_polydesc.pverts;

	for (i=0 ; i<r_polydesc.nump ; i++)
	{
		if (pverts->v < ymin)
		{
			ymin = pverts->v;
			s_minindex = i;
		}

		if (pverts->v > ymax)
		{
			ymax = pverts->v;
			s_maxindex = i;
		}

		pverts++;
	}

	ymin = ceil (ymin);
	ymax = ceil (ymax);

	if (ymin >= ymax)
		return; // doesn't cross any scans at all

	cachewidth = r_polydesc.pixel_width;
	cacheheight = r_polydesc.pixel_height;
	cacheblock = r_polydesc.pixels;

	// copy the first vertex to the last vertex, so we don't have to deal with
	// wrapping
	nump = r_polydesc.nump;
	pverts = r_polydesc.pverts;
	pverts[nump] = pverts[0];

	// PS1 style - calculate affine coefficients
	R_PolygonCalculateAffineCoefficients(pverts);
	
	// Calculate Z gradients only (for depth buffer)
	R_PolygonCalculateGradients(&d_ziorigin, &d_zistepu, &d_zistepv);
	
	// Scan edges and draw
	R_PolygonScanLeftEdge(spans);
	R_PolygonScanRightEdge(spans);

	R_PolygonDrawSpans(spans, iswater, d_ziorigin, d_zistepu, d_zistepv);
}

/*
** R_DrawAlphaSurfaces
*/
void
R_DrawAlphaSurfaces(const entity_t *currententity)
{
	msurface_t *s = r_alpha_surfaces;
	const model_t *currentmodel = r_worldmodel;

	modelorg[0] = -r_origin[0];
	modelorg[1] = -r_origin[1];
	modelorg[2] = -r_origin[2];

	while ( s )
	{
		R_BuildPolygonFromSurface(currententity, currentmodel, s);

		// pass down all the texinfo flags, not just SURF_WARP.
		if (s->texinfo->flags & SURF_TRANS66)
			R_ClipAndDrawPoly( 0.60f, (s->texinfo->flags & (SURF_WARP | SURF_SCROLL)), true );
		else
			R_ClipAndDrawPoly( 0.30f, (s->texinfo->flags & (SURF_WARP | SURF_SCROLL)), true );

		s = s->nextalphasurface;
	}

	r_alpha_surfaces = NULL;
}

/*
** R_IMFlatShadedQuad
*/
void
R_IMFlatShadedQuad( const vec3_t a, const vec3_t b, const vec3_t c, const vec3_t d, int color, float alpha )
{
	vec3_t s0, s1;

	r_polydesc.nump = 4;
	VectorCopy( r_origin, r_polydesc.viewer_position );

	VectorCopy( a, r_clip_verts[0][0] );
	VectorCopy( b, r_clip_verts[0][1] );
	VectorCopy( c, r_clip_verts[0][2] );
	VectorCopy( d, r_clip_verts[0][3] );

	r_clip_verts[0][0][3] = 0;
	r_clip_verts[0][1][3] = 0;
	r_clip_verts[0][2][3] = 0;
	r_clip_verts[0][3][3] = 0;

	r_clip_verts[0][0][4] = 0;
	r_clip_verts[0][1][4] = 0;
	r_clip_verts[0][2][4] = 0;
	r_clip_verts[0][3][4] = 0;

	VectorSubtract( d, c, s0 );
	VectorSubtract( c, b, s1 );
	CrossProduct( s0, s1, r_polydesc.vpn );
	VectorNormalize( r_polydesc.vpn );

	r_polydesc.dist = DotProduct( r_polydesc.vpn, r_clip_verts[0][0] );

	r_polyblendcolor = color;

	R_ClipAndDrawPoly( alpha, false, false );
}
