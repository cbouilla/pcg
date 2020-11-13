from itertools import islice
import random

# actual implementation of the mystery function.
state = 0x123456deadbeefbadcafe
for i in range(16):
    state = (state * 0xda942042e4dd58b5) % (2**128)
    print("{}   :   {:016x}".format(state, state >> 64))
    if i == 7:
        print('----------------------------')
    
print('=================================')


def mystery(a, b, c, d, e, f, g, h):
    """
    >>> seed = [0xa872551a0848efe8, 0xeb568cb9eddd0de1, 0x72582f39ab1530c7, 0xf13af01e3a960ad8, \
                0xc6bfb791661c6d3f, 0x71c294663911d20e, 0x41dc534da47e7130, 0x141e26c46484e8f2]
    >>> stream = mystere(*seed)
    >>> for x in islice(stream, 0, 8):
    ...     print(hex(x))
    0x3f76377c69c0d314
    0xc8a788c2474ec97
    0x8ee93c62db44939b
    0xf385b6d382ab0e42
    0x597e35963c650290
    0xd8ca9a66f7b08d1a
    0x2b05a43500bcbe58
    0x157cd33ef08dfeac
    """
    prime = 2**64 - 59
    while True:
        i =   7886*a  +7612*b +26661*c -15282*d +32867*e +48473*f -36291*g +35350*h
        j = -57569*a +19465*b -24992*c  -1812*d +26618*e +29956*f +34577*g -18859*h
        k =  -2446*a +61735*b -22427*c +34186*d  +6992*e  -9478*f  +8197*g -27775*h
        l =  -1185*a  -6312*b -56926*c +15565*d +11059*e  +1553*f -57153*g  +5707*h
        m =   -535*a  -1791*b -39310*c +43305*d  +6273*e +56313*f  -3904*g +17665*h
        n = -10209*a +22995*b +20231*c  +1049*d -14251*e -41664*f +12854*g +33455*h
        o =  23582*a -29517*b -16622*c +16105*d +39186*e +13879*f +16248*g +11659*h
        p =  11671*a  +5093*b -13453*c -33945*d  +6320*e +22666*f +28441*g +25619*h
        q = round(i / prime)
        r = round(j / prime)
        s = round(k / prime)
        t = round(l / prime)
        u = round(m / prime)
        v = round(n / prime)
        w = round(o / prime)
        x = round(p / prime)
        y = 0
        y = y -  100539497500637377139787102742151 * q
        y = y - 4385199133127559891011963479249808 * r
        y = y + 2714301933393960341823034119028324 * s
        y = y -  893200328965818184420614031056548 * t
        y = y -  392009643508364887148804587723243 * u
        y = y - 2421047828928147225283878252264788 * v
        y = y + 1566125033113415705045963357974976 * w
        y = y + 2771473709774666053247393143332435 * x
        y = y * 78835592693290935975773860076245338465
        z = (y >> 64) & 0xffffffffffffffff
        a = b
        b = c
        c = d
        d = e
        e = f
        f = g
        g = h
        h = z
        yield z
