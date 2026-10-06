from math import ulp


class ConvexPolygon:
    """Convex polygon in the plane, built from a list of vertices (x, y).

    Vertices may be given clockwise or counterclockwise. Three neighbouring
    vertices on one straight line are rejected.
    """

    def __init__(self, verts):
        if not ConvexPolygon.is_valid_format(verts):
            raise ValueError("need a list of at least 3 vertices, each (x, y) with numbers")
        if not ConvexPolygon.is_convex(verts):
            raise ValueError("vertices do not form a convex polygon")
        self._verts = [tuple(vert) for vert in verts]

    @property
    def vertices(self):
        return list(self._verts)

    def __repr__(self):
        return "ConvexPolygon(" + str(self._verts) + ")"

    #methods

    def perimeter(self):
        n = len(self._verts)
        total = 0
        for i in range(n):
            total += ConvexPolygon.distance(self._verts[i], self._verts[(i + 1) % n])
        return total

    def area(self):
        n = len(self._verts)
        total = 0
        for i in range(n):
            total += ConvexPolygon.cross(self._verts[i], self._verts[(i + 1) % n])
        return abs(total) / 2

    #operations

    def contains(self, other):
        """True if a point (x, y) or a whole ConvexPolygon lies inside (boundary counts)."""
        if isinstance(other, ConvexPolygon):
            for p in other.vertices:
                if not self.contains(p):
                    return False
            return True
        n = len(self._verts)
        orientation = ConvexPolygon.turn(self._verts[0], self._verts[1], self._verts[2])
        for i in range(n):
            s = ConvexPolygon.side(self._verts[i], self._verts[(i + 1) % n], other)
            if s != 0 and (s > 0) != (orientation > 0):
                return False
        return True

    def triangulate(self):
        """Split into n - 2 triangles fanning out from the first vertex."""
        triangles = []
        for i in range(1, len(self._verts) - 1):
            triangles.append(ConvexPolygon([self._verts[0], self._verts[i], self._verts[i + 1]]))
        return triangles

    def intersection(self, other):
        """Shared region of two convex polygons, or None if they share no area."""
        inside_sign = 1 if ConvexPolygon.turn(other._verts[0], other._verts[1], other._verts[2]) > 0 else -1
        result = list(self._verts)
        m = len(other._verts)
        for j in range(m):
            if not result:
                break
            result = ConvexPolygon.clip(result, other._verts[j], other._verts[(j + 1) % m], inside_sign)
        result = ConvexPolygon.clean(result)
        if not result:
            return None
        return ConvexPolygon(result)

    #static helpers 

    @staticmethod
    def is_valid_format(verts):
        if not isinstance(verts, (list, tuple)):
            return False
        if len(verts) < 3:
            return False
        for vert in verts:
            if not isinstance(vert, (list, tuple)):
                return False
            if len(vert) != 2:
                return False
            if not (isinstance(vert[0], (int, float)) and isinstance(vert[1], (int, float))):
                return False
        return True

    @staticmethod
    def cross(u, v):
        return u[0] * v[1] - u[1] * v[0]

    @staticmethod
    def turn(a, b, c):
        # > 0: left turn at b, < 0: right turn, 0: straight
        return ConvexPolygon.side(a, b, c)

    @staticmethod
    def side(a, b, p):
        # > 0: p is left of the line a -> b, < 0: right, 0: on the line
        left = (b[0] - a[0]) * (p[1] - a[1])
        right = (b[1] - a[1]) * (p[0] - a[0])
        result = left - right
        # Allow a few floating-point rounding steps when the products cancel.
        # This follows their scale, so small polygons are not rounded away.
        if isinstance(result, float) and abs(result) <= 4 * (ulp(left) + ulp(right)):
            return 0
        return result

    @staticmethod
    def distance(p, q):
        dx = q[0] - p[0]
        dy = q[1] - p[1]
        return (dx ** 2 + dy ** 2) ** 0.5

    @staticmethod
    def is_convex(verts):
        n = len(verts)
        first_turn = ConvexPolygon.turn(verts[0], verts[1], verts[2])
        if first_turn == 0:
            return False
        for i in range(n):
            a = verts[i]
            b = verts[(i + 1) % n]
            # Every other vertex must be strictly on the inner side of this edge.
            # Checking only neighbouring turns would also accept some stars.
            for j in range(n):
                if j == i or j == (i + 1) % n:
                    continue
                s = ConvexPolygon.side(a, b, verts[j])
                if s == 0 or (s > 0) != (first_turn > 0):
                    return False
        return True

    @staticmethod
    def line_crossing(p, q, a, b):
        # Point where segment p -> q crosses the line through a and b
        p_side = ConvexPolygon.side(a, b, p)
        q_side = ConvexPolygon.side(a, b, q)
        if p_side == 0:
            return tuple(p)
        if q_side == 0:
            return tuple(q)
        t = p_side / (p_side - q_side)
        return (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))

    @staticmethod
    def clip(subject, a, b, inside_sign):
        # Keep the part of `subject` on the inner side of the line a -> b (Sutherland-Hodgman step)
        result = []
        n = len(subject)
        for i in range(n):
            p = subject[i]
            q = subject[(i + 1) % n]
            p_in = ConvexPolygon.side(a, b, p) * inside_sign >= 0
            q_in = ConvexPolygon.side(a, b, q) * inside_sign >= 0
            if p_in:
                result.append(p)
            if p_in != q_in:
                result.append(ConvexPolygon.line_crossing(p, q, a, b))
        return result

    @staticmethod
    def clean(verts):
        # Drop repeated points and straight vertices without rounding coordinates.
        points = []
        for p in verts:
            p = tuple(p)
            if p not in points:
                points.append(p)
        changed = True
        while changed and len(points) >= 3:
            changed = False
            n = len(points)
            for i in range(n):
                if ConvexPolygon.turn(points[i - 1], points[i], points[(i + 1) % n]) == 0:
                    points.pop(i)
                    changed = True
                    break
        if len(points) < 3:
            return []
        return points


square = ConvexPolygon([(0, 0), (4, 0), (4, 4), (0, 4)])
shifted = ConvexPolygon([(2, 2), (6, 2), (6, 6), (2, 6)])
small = ConvexPolygon([(1, 1), (3, 1), (3, 3), (1, 3)])
far = ConvexPolygon([(10, 10), (12, 10), (12, 12)])

print(square.vertices)                 # [(0, 0), (4, 0), (4, 4), (0, 4)]
print(square.area())                   # 16.0
print(square.perimeter())              # 16.0

print(square.contains((2, 2)))         # True
print(square.contains((4, 2)))         # True  (on the edge)
print(square.contains((5, 2)))         # False
print(square.contains(small))          # True
print(square.contains(shifted))        # False

print(square.triangulate())            # 2 triangles

common = square.intersection(shifted)
print(common, common.area())           # 2x2 square, 4.0
print(square.intersection(far))        # None

for bad in [[(0, 0), (4, 0), (2, 2), (4, 3), (0, 3)],   # dent
            [(0, 0), (2, 0), (4, 0), (4, 3), (0, 3)],   # straight vertex
            [(0, 0), (4, 0)],                           # too few vertices
            [(0, 0), (4, "x"), (0, 3)]]:                # bad format
    try:
        ConvexPolygon(bad)
    except ValueError as error:
        print("Rejected:", error)

try:
    square.vertices = [(0, 0), (1, 0), (0, 1)]
except AttributeError:
    print("vertices is read-only")
