#ifndef __EVENTS_H_
#define __EVENTS_H_

#include "triangle.h"

/**
 * \brief Scattering event in the bulk of a material
 */
struct scatter_event
{
	/**
	 * \brief Type of event
	 *
	 * 0 signifies no event; otherwise, an index to the event type.
	 */
	uint8_t type;

	/**
	 * \brief Distance to the event.
	 */
	real distance;
};

/**
 * \brief Intersection event with geometry
 */
struct intersect_event
{
	/**
	 * \brief Distance to the event.
	 */
	real isect_distance;

	/**
	 * \brief Pointer to the triangle intersected.
	 */
	triangle* isect_triangle;
};

#endif // __EVENTS_H_
