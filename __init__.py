# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 3
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

bl_info = {
    "name": "Ray Caster",
    "description":  "Use the default light as a ray caster",
    "author": "Noizirom",
    "version": (0, 0, 1),
    "blender": (2, 83, 0),
    "location": "View3D > UI > Ray Caster",
    "warning": "", 
    "wiki_url": "",
    "tracker_url": "",
    "category": "Ray Caster"
}

import bpy


def dome_ops(rad, seg, rings):
    """
    Calculates a dome object geometry from 
    args:
        rad: radius length
        seg: segment count
        rings: ring count
    return:
        verts: array of vertices
        edges: array of edges
        faces: array of faces vertices indexes
    """
    from numpy import array, empty, vstack, radians, cos, sin
    count = seg * rings
    iseg = seg-1
    icount = count-1
    last_ring = count-seg
    verts, edges = [], []
    s_ang = 360 / seg
    r_ang = 90 / rings
    idx = [i for i in range(seg) if i % 2 == 0]
    for r in range(rings):
        r_theta = r * r_ang
        r_theta = radians(r_theta)
        r_rad = cos(r_theta) * rad
        step = r * seg
        e_list = [[i + step, (i + 1) + step] for i in range((seg - 1))]
        e_end = [e_list[-1][1], e_list[0][0]]
        e_copy = e_list.copy()
        e_copy.append(e_end.copy())
        if r != 0 and r % 2 != 0:
            e_list2 = e_list[::-1].copy()
            e_list2 = [[i[1], i[0]] for i in e_list2]
            e_end = [e_list2[-1][1], e_list2[0][0]]
            edges.append(e_end)
            edges = edges + e_list2
        else:
            edges = edges + e_list
            edges.append(e_end)
        if r > 0:
            sides = []
            for i in e_copy:
                side = [[i[0], i[0] - seg], [i[1] - seg, i[1]]]
                sides.append(side)
            for s in idx:
                edges = edges + sides[s]
        if r == (rings - 1):
            e_top =[[[count, i[0]], [i[1], count]] for i in e_copy]
            i_list = []
            for i in idx:
                i_list = i_list + e_top[i]
            edges = edges + i_list
        for s in range(seg):
            s_theta = s * s_ang
            s_theta = radians(s_theta)
            z = float(round(((sin(r_theta)) * rad), ndigits=5))
            co = [float(round((cos(s_theta) * (rad*cos(r_theta))), ndigits=5)),
                float(round((sin(s_theta) * (rad*cos(r_theta))), ndigits=5)), z]
            verts.append(co)
    if rings > 1:
        verts.append([0.,0., rad])
    else:
        verts.append([0.,0.,0.])
    f_close = []
    for r in range(rings-1):
        segr = seg*r
        f_close.append([iseg+segr, 0+segr, seg+segr, seg+iseg+segr])
    f_top_idx = [i for i in range((last_ring),count)]
    ft = [[i, i+1, count] for i in f_top_idx if i != (icount)]
    f_top = ft + [[icount,last_ring,count]]
    f_sides = [[i, i + 1, i + 1 + seg, i + seg] for i in range(last_ring) if i not in [(r*seg-1) for r in range(rings) if r!=0]]
    faces = f_sides + f_close + f_top
    verts = array(verts)
    vct = len(verts)
    vts = empty(vct*3).reshape((vct, 3))
    vts[:,0] = -verts[:,1]
    vts[:,1] = verts[:,0]
    vts[:,2] = verts[:,2]
    return vts.tolist(), edges, faces

def search_orb(target_ob, origin, segment_count, depsgraph=None, seperate=False):
    """
    Creates a spherical ray casting object 
    args:
        target_ob:     object to cast rays at
        origin:        origin point of ray cast
        segment_count: a factor used to set 
                       segment and ring count for ray vectors
                       total number of rays = (segment_count * segment_count + 1) * 2 - segment_count
        depsgraph:     link to depsgraph
        seperate:      return raycast_pos and raycast_neg as seperate arrays if True, A combined array if False
    return:
        raycast_pos: array of results of ray casts in positive direction from the origin xy plane
                success:   Bool of target intersection
                hit_point: Point of intersection if success is True, [0,0,0] if success is False
                hit_norm:  Normal vector of intersection if success is True, [0,0,0] if success is False
                hit_index: Index of polygon intersection if success is True, -1 if success is False
        raycast_neg: array of results of ray casts in negative direction from the origin xy plane
                success:   Bool of target intersection
                hit_point: Point of intersection if success is True, [0,0,0] if success is False
                hit_norm:  Normal vector of intersection if success is True, [0,0,0] if success is False
                hit_index: Index of polygon intersection if success is True, -1 if success is False
    """
    from numpy import array, vstack
    def ray_sphere(segment_count):
        rays = array(dome_ops(1, segment_count, segment_count)[0])
        negative_rays = -rays[segment_count:]
        return rays, negative_rays
    def results(inspect):
        result = []
        for ray in inspect:
            success, hit_point, hit_normal, hit_index = cast(ray)
            result.append([success, (array(hit_point) if success else ray), array(hit_normal), hit_index])
        return array(result)
    ray_inspect, ray_inspect_neg = ray_sphere(segment_count)
    cast = lambda direction: target_ob.ray_cast(origin, direction, depsgraph=depsgraph)
    raycast_pos = results(ray_inspect)
    raycast_neg = results(ray_inspect_neg)
    if seperate:
        return raycast_pos, raycast_neg
    else:
        return vstack((raycast_pos, raycast_neg))

def hit_indexes(search):
    from numpy import argwhere
    return argwhere(search[:,0])[:,0]

def miss_indexes(search):
    from numpy import argwhere, invert
    return argwhere(search[:,0] == False)[:,0]

class DrawRays:
    def __init__(self, coords, orb_loc, color=(1,1,1,1)):
        self.coords = coords
        self.count = len(self.coords)
        self.orb_loc = orb_loc
        self.color = color
        self.handle = None
    #
    def create(self):
        import gpu
        from gpu_extras.batch import batch_for_shader
        from numpy import array, arange, empty, zeros
        rays = array([self.orb_loc] + self.coords.tolist())
        r_ct = rays.shape[0]
        indices = empty((r_ct-1, 2), dtype=int)
        indices[:,0] = zeros(r_ct-1, dtype=int)
        indices[:,1] = arange(r_ct)[1:]
        indices = indices.tolist()
        shader = gpu.shader.from_builtin('3D_UNIFORM_COLOR')
        batch = batch_for_shader(shader, 'LINES', {"pos": self.coords.tolist()}, indices=indices)
        shader.bind()
        shader.uniform_float("color", self.color)
        batch.draw(shader)
        self.update()
    #
    def draw(self):
        if self.handle != None:
            self.remove()
        self.handle = bpy.types.SpaceView3D.draw_handler_add(self.create, (), "WINDOW", 'POST_VIEW')
    #
    def remove(self):
        bpy.types.SpaceView3D.draw_handler_remove(self.handle, 'WINDOW')
        self.handle = None
    #
    def update(self):
        for area in bpy.context.window.screen.areas:
            if area.type == 'VIEW_3D':
                area.tag_redraw()

class RayCasterOrb:
    def __init__(self, target_ob, orb, ray_factor=2, HIT_COLOR=(0,1,0,1), MISS_COLOR=(1,0,0,1), CONNECT_COLOR=(0,0,1,1)):
        self.target_ob = target_ob
        self.orb = orb
        self.ray_factor = ray_factor
        self.HIT_COLOR = HIT_COLOR
        self.MISS_COLOR = MISS_COLOR
        self.CONNECT_COLOR = CONNECT_COLOR
        self.hit_rays = None
        self.miss_rays = None
        self.connection_ray = None
    #
    def draw(self):
        from numpy import array
        orb_loc = array(self.orb.location)
        targ_loc = array(self.target_ob.location)
        search = search_orb(self.target_ob, orb_loc, self.ray_factor, depsgraph=bpy.context.evaluated_depsgraph_get(), seperate=False)
        hit = hit_indexes(search)
        miss = miss_indexes(search)
        self.connection_ray = DrawRays(array([targ_loc]), orb_loc, self.CONNECT_COLOR)
        self.connection_ray.draw()
        coords = array([orb_loc] + [(i-targ_loc).tolist() for i in search[hit][:,1]])
        mcoords = array([orb_loc] + [(i+orb_loc).tolist() for i in search[miss][:,1] if i.tolist()])
        if hit.size > 0 and bpy.context.scene.ray_props.show_hits:
            #print(f"Hit points: {coords[1:]}")
            self.hit_rays = DrawRays(coords, orb_loc, self.HIT_COLOR)
            self.hit_rays.draw()
        if miss.size > 0 and bpy.context.scene.ray_props.show_misses:
            #print(f"Missed Rays: {coords}")
            self.miss_rays = DrawRays(mcoords, orb_loc, self.MISS_COLOR)
            self.miss_rays.draw()
    #
    def remove(self):
        try:
            self.connection_ray.remove()
            self.connection_ray = None
        except Exception as e:
                print(e)
        try:
            if self.hit_rays != None:
                self.hit_rays.remove()
        except Exception as e:
                print(e)
        try:
            if self.miss_rays != None:
                self.miss_rays.remove()
        except Exception as e:
                print(e)
    #
    def update(self):
        if self.connection_ray != None:
            try:
                self.remove()
            except Exception as e:
                print(e)
        self.draw()


#####################################################################

class RayProp(bpy.types.PropertyGroup):
    #
    show_hits: bpy.props.BoolProperty(
        name = "Show_Hit_Rays",
        description = "Show rays where hit",
        default = True,
        )
    #
    show_misses: bpy.props.BoolProperty(
        name = "Show_Miss_Rays",
        description = "Show rays where missed",
        default = True,
        )
    #
    ray_fac: bpy.props.IntProperty(
        name="Ray_Cast_factor",
        default=2,
        soft_min=2,
        soft_max=100,
        description="Factor for determining amount of segments and ring points that the rays use")
    #
    target_o: bpy.props.StringProperty(
        name="Target Object",
        default="Cube",
        description="Target object to cast rays on")

raycasterorb = None

class ModalTimerOperator(bpy.types.Operator):
    """Operator which runs its self from a timer"""
    bl_idname = "wm.modal_timer_operator"
    bl_label = "Modal Timer Operator"

    _timer = None

    def modal(self, context, event):
        global raycasterorb
        if event.type in {'RIGHTMOUSE', 'ESC'}:
            raycasterorb.remove()
            self.cancel(context)
            return {'CANCELLED'}
        if event.type == 'TIMER':
            raycasterorb.update()
        return {'PASS_THROUGH'}

    def execute(self, context):
        wm = context.window_manager
        self._timer = wm.event_timer_add(0.1, window=context.window)
        wm.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def cancel(self, context):
        wm = context.window_manager
        raycasterorb.remove()
        wm.event_timer_remove(self._timer)


class CastRays(bpy.types.Operator):
    """Cast rays from a single point"""
    bl_idname = "rayc.ray_caster"
    bl_label = "Ray Caster"
    bl_options = {'REGISTER', 'INTERNAL', 'UNDO'}
    
    @classmethod
    def poll(cls, context):
        return context.mode in {'OBJECT'} and context.object == bpy.data.objects['Light']
    
    def execute(self, context):
        global raycasterorb
        scn = bpy.context.scene
        ray_factor=scn.ray_props.ray_fac
        target_object = scn.ray_props.target_o
        raycasterorb = RayCasterOrb(bpy.data.objects.get(target_object, None), bpy.data.objects["Light"], ray_factor=ray_factor, HIT_COLOR=(0,1,0,1), MISS_COLOR=(1,0,0,1), CONNECT_COLOR=(0,0,1,1))
        if scn.ray_props.show_hits:
            bpy.ops.wm.modal_timer_operator()
        return {'FINISHED'}

def draw_rays_layout(self, context, layout):
    scn = context.scene
    ray_box = layout.box()
    ray_box.label(text="RayCaster")
    ray_box.operator(CastRays.bl_idname, icon="USER")
    if context.object == bpy.data.objects["Light"]:
        ray_box.prop(scn.ray_props ,"target_o")
        ray_box.prop(scn.ray_props ,"ray_fac")
        ray_box.prop(scn.ray_props ,"show_hits")
        ray_box.prop(scn.ray_props ,"show_misses")


class OBJECT_PT_RayCPanel(bpy.types.Panel):
    bl_label = f"Ray Caster"
    bl_idname = "OBJECT_PT_raycast_panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Ray Caster"
    bl_context = "objectmode"
    #
    def draw(self, context):
        layout = self.layout
        draw_rays_layout(self, context, layout)


#Registry
####################

classes = [
    RayProp,
    ModalTimerOperator,
    CastRays,
    OBJECT_PT_RayCPanel,
]

rev_classes = classes.copy()
rev_classes.reverse()

def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.Scene.ray_props = bpy.props.PointerProperty(type=RayProp)


def unregister():
    for cls in rev_classes:
        bpy.utils.unregister_class(cls)
    #
    del bpy.types.Scene.ray_props


if __name__ == "__main__":
    register()

