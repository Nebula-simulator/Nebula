/**
 * @file src/cdsem/tribox.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__TRIBOX__HEADER_INCLUDED
#define eSCATTER__CDSEM__TRIBOX__HEADER_INCLUDED

int triBoxOverlap(double boxcenter[3], double boxhalfsize[3], double triverts[3][3]);

#endif
