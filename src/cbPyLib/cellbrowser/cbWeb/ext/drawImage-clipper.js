/**
 * MIT License - Copyright (c) 2021 Kaiido
 *
 * A monkey-patch for Safari's drawImage.
 *
 * This browser doesn't handle well using the cropping abilities of drawImage
 * with out-of-bounds values.
 * (see https://stackoverflow.com/questions/35500999/cropping-with-drawimage-not-working-in-safari)
 * This script takes care of detecting when the monkey-patch is needed,
 * and does redefine the cropping parameters so they fall inside the source's boundaries.
 *
**/

(()=> {

  if( !needPoly() ) { return; }

  const proto = CanvasRenderingContext2D.prototype;
  const original = proto.drawImage;
  if( !original ) {
    console.error( "This script requires a basic implementation of drawImage" );
    return;
  }

  proto.drawImage = function drawImage( source, x, y ) { // length: 3

    const will_crop = arguments.length === 9;
    if( !will_crop ) {
      return original.apply( this, [...arguments] );
    }

    const safe_rect = getSafeRect( ...arguments );
    if( isEmptyRect( safe_rect ) ) {
      return;
    }
    return original.apply( this, safe_rect );
  } 

  function needPoly() {
    const ctx = document.createElement( "canvas" ).getContext( "2d" );
    ctx.fillRect( 0, 0, 40, 40 );
    ctx.drawImage( ctx.canvas, -40, -40, 80, 80, 50, 50, 20, 20 );

    const img = ctx.getImageData( 50, 50, 30, 30 ); // 10px around expected square
    const data = new Uint32Array( img.data.buffer );
    const colorAt = (x, y) => data[ y * img.width + x ];

    const transparents = [ [ 9, 9 ], [ 20, 9 ], [ 9, 20 ], [ 20, 20 ] ];
    const blacks = [ [ 10, 10 ], [ 19, 10 ], [ 10, 19 ], [ 19, 19 ] ];
    return transparents.some( ([ x, y ]) => colorAt( x, y ) !== 0x00000000 ) ||
      blacks.some( ([ x, y ]) => colorAt( x, y ) === 0x00000000 )
  }

  function getSafeRect( image, sx, sy, sw, sh, dx, dy, dw, dh ) {
  
    const { width, height } = getSourceDimensions( image );
    
    if( sw < 0 ) {
      sx += sw;
      sw = Math.abs( sw );
    }
    if( sh < 0 ) {
      sy += sh;
      sh = Math.abs( sh );
    }
    if( dw < 0 ) {
      dx += dw;
      dw = Math.abs( dw );
    }
    if( dh < 0 ) {
      dy += dh;
      dh = Math.abs( dh );
    }
    const x1 = Math.max( sx, 0 );
    const x2 = Math.min( sx + sw, width );
    const y1 = Math.max( sy, 0 );
    const y2 = Math.min( sy + sh, height );
    const w_ratio = dw / sw;
    const h_ratio = dh / sh;

    return [
      image,
      x1,
      y1,
      x2 - x1,
      y2 - y1,
      sx < 0 ? dx - (sx * w_ratio) : dx,
      sy < 0 ? dy - (sy * h_ratio) : dy,
      (x2 - x1) * w_ratio,
      (y2 - y1) * h_ratio
    ];

  }

  function isEmptyRect( args ) {
    // sw, sh, dw, dh
    return [ 3, 4, 7, 8 ].some( (index) => !args[ index ] );
  }

  function getSourceDimensions( source ) {
    const sourceIs = ( type ) => {
      const constructor = globalThis[ type ];
      return constructor && (source instanceof constructor);
    };
    if( sourceIs( "HTMLImageElement" ) ) {
      return { width: source.naturalWidth, height: source.naturalHeight };
    }
    else if( sourceIs( "HTMLVideoElement" ) ) {
      return { width: source.videoWidth, height: source.videoHeight };
    }
    else if( sourceIs( "SVGImageElement" ) ) {
      throw new TypeError( "SVGImageElement isn't yet supported as source image.", "UnsupportedError" );
    }
    else if( sourceIs( "HTMLCanvasElement" ) || sourceIs( "ImageBitmap" ) ) {
      return source;
    }
  }

})();
