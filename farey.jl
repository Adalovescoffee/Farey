### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 42a780fe-1659-4d9b-b900-3321ce025dd6
using Luxor, PlutoUI, Colors, SimpleWeightedGraphs, Graphs, GraphPlot

# ╔═╡ b03a4559-a960-4b56-97ed-a144407ba2eb
using GeometryTypes

# ╔═╡ 15f0e8a0-59e4-11ef-288b-17fd5d6af339
md" # Implementing farey diagrams in julia  
I am going to be using graphs.jl and luxor package 
"

# ╔═╡ 0dfb78a9-a549-4c56-92f5-f4e5aa6f26d5
md"
Setting up an object that has (,) as values so it's going to be my points 
"

# ╔═╡ 68057bda-5823-4e28-8a56-cfdfedb15e5d
struct FareyVertex 
		value::Set{Rational{Int}}
end 

# ╔═╡ 70d7cb70-f755-4241-bbbd-77852f0b99d5
function create_custom_graph()
	g = SimpleWeightedGraph(0)
	vertex_value = Dict{Int, Set{Rational{Int}}}()
	return g, vertex_value
end

# ╔═╡ 817a2afb-3845-4cc7-bf9d-cb6f78b50de6
function add_farey_vertex!(g::SimpleWeightedGraph, vertex_value::Dict{Int, Set{Rational{Int}}}, values::Set{Rational{Int}})
	#nv is for number of vertices in ur initial graph
    v_index = nv(g) + 1
    add_vertices!(g, 1)  # Add one vertex to the graph
    vertex_value[v_index] = values
    return v_index
end


# ╔═╡ ded767d6-493f-40af-bba6-7f3bd5f9b57a
md"
Add_farey_vertex takes a graph g, its vertexvalues which is a dictionnary
"

# ╔═╡ 2c41b5bc-794a-432e-8076-5be96de4bc70
function add_farey_edge!(g::SimpleWeightedGraph, v1::Int, v2::Int,weight::Float64)
    add_edge!(g, v1,v2,weight)
end
# I suppose i'll get the v1 and v2 from the previous point feel like this is useless we'll see idk :( 

# ╔═╡ b0985351-2382-4ecb-9bb1-1ae8bb9271a9
function update_edge_weight!(g::SimpleWeightedGraph, v1::Int, v2::Int, new_weight::Float64)
    # This function updates the weight of an existing edge
    if has_edge(g, v1, v2)
        g.weights[v1, v2] = new_weight
        g.weights[v2, v1] = new_weight  # Ensure symmetry
    end
end

# ╔═╡ a6989575-1ced-42ec-bb25-dcde4a750b71
md"## Initializing Farey diagram for n = 1 "

# ╔═╡ 8183a1a1-af93-4e48-bf57-7d1443a9e79e
function FareyGraph(n)
    if n == 1
        F_1, Vval = create_custom_graph()
        v1 = add_farey_vertex!(F_1, Vval, Set(Rational(1, 1)))
        v2 = add_farey_vertex!(F_1, Vval, Set([Rational(1, 0), Rational(-1, 0)]))  # Handles ∞ case
        v3 = add_farey_vertex!(F_1, Vval, Set(Rational(-1, 1)))
        v4 = add_farey_vertex!(F_1, Vval, Set(Rational(0, 1)))
        
        add_farey_edge!(F_1, v1, v2, 0.1)
        add_farey_edge!(F_1, v2, v3, 0.1)
        add_farey_edge!(F_1, v3, v4, 0.1)
        add_farey_edge!(F_1, v1, v4, 0.1) 
        return F_1, Vval
        
    else 
        F = FareyGraph(n-1)
        Graph_n, Value_n = F
        CopyG = SimpleWeightedGraph(Graph_n)  # Create a new graph instead of copying
        for e in edges(CopyG)
            p1, p2 = src(e), dst(e)
            w = weight(e)
            f_1 = first(Value_n[p1])
            f_2 = first(Value_n[p2])
            if f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
                if f_1 in [Rational(1, 0), Rational(-1, 0)]
                	if f_1 == Rational(1,0)
						f_2<0 ? f_1 = -f_1 : f_1 = f_1
					else 
						f_2>0 ? f_1 = -f_1 : f_1 = f_1
					end
				else
                    if f_2 == Rational(1,0)
						f_1<0 ? f_2 = -f_2 : f_2 = f_2
					else 
						f_1>0 ? f_2 = -f_2 : f_2 = f_2
					end
                end
            end
            new_value = Rational(numerator(f_1) + numerator(f_2), denominator(f_1) + denominator(f_2))
            if w ≈ 0.1
                v_child = add_farey_vertex!(Graph_n, Value_n, Set(new_value))
                add_farey_edge!(Graph_n, v_child, p1, 0.1)  
                add_farey_edge!(Graph_n, v_child, p2, 0.1)
                update_edge_weight!(Graph_n, p1, p2, 1.0)  # Remove the old edge instead of setting weight to 0
            end 
        end
        return Graph_n, Value_n
    end
end

# ╔═╡ d5de3fdc-ec9d-40ac-8ba2-ff1db06e8c72
function FareyGraph2(n)
    if n == 1
        F_1, Vval = create_custom_graph()
        v1 = add_farey_vertex!(F_1, Vval, Set(Rational(1, 1)))
        v2 = add_farey_vertex!(F_1, Vval, Set([Rational(1, 0), Rational(-1, 0)]))  # Handles ∞ case
        v3 = add_farey_vertex!(F_1, Vval, Set(Rational(-1, 1)))
        v4 = add_farey_vertex!(F_1, Vval, Set(Rational(0, 1)))
        v5 = add_farey_vertex!(F_1, Vval, Set(Rational(2/1)))
		v6 = add_farey_vertex!(F_1, Vval, Set(Rational(-2/1)))
		v7 = add_farey_vertex!(F_1,Vval, Set(Rational(1/2)))
		v8 = add_farey_vertex!(F_1,Vval,Set(Rational(-1/2)))
        add_farey_edge!(F_1, v1, v3, 0.1)
        add_farey_edge!(F_1, v2, v7, 0.1)
		add_farey_edge!(F_1,v4,v5, 0.1)
		add_farey_edge!(F_1,v4,v6, 0.1)
		add_farey_edge!(F_1,v2,v8, 0.1)
        return F_1, Vval
        
    else 
        F = FareyGraph2(n-1)
        Graph_n, Value_n = F
        CopyG = SimpleWeightedGraph(Graph_n)  # Create a new graph instead of copying
        for e in edges(CopyG)
            p1, p2 = src(e), dst(e)
            w = weight(e)
            f_1 = first(Value_n[p1])
            f_2 = first(Value_n[p2])
			
				
			
		if f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
                if f_1 in [Rational(1, 0), Rational(-1, 0)]
                	if f_1 == Rational(1,0)
						f_2<0 ? f_1 = -f_1 : f_1 = f_1
					else 
						f_2>0 ? f_1 = -f_1 : f_1 = f_1
					end
				else
                    if f_2 == Rational(1,0)
						f_1<0 ? f_2 = -f_2 : f_2 = f_2
					else 
						f_1>0 ? f_2 = -f_2 : f_2 = f_2
					end
                end
		
		end
            
			
            if w ≈ 0.1
				new_value1 = Rational(numerator(f_1) + 2*numerator(f_2), denominator(f_1) + 2*denominator(f_2))
				new_value2 = Rational(2*numerator(f_1) + numerator(f_2), 2*denominator(f_1)  + denominator(f_2))
                v_child1 = add_farey_vertex!(Graph_n, Value_n, Set(new_value1))
				v_child2 = add_farey_vertex!(Graph_n, Value_n, Set(new_value2))
				new_value3 = Rational(numerator(f_1) - 2*numerator(f_2), denominator(f_1) + -2*denominator(f_2))
				new_value4 = Rational(-2*numerator(f_1) + numerator(f_2), -2*denominator(f_1)  + denominator(f_2))
                v_child3 = add_farey_vertex!(Graph_n, Value_n, Set(new_value3))
				v_child4 = add_farey_vertex!(Graph_n, Value_n, Set(new_value4))
				add_farey_edge!(Graph_n, v_child2, p1, 0.1)  
                add_farey_edge!(Graph_n, v_child1, p2, 0.1)
                add_farey_edge!(Graph_n, v_child4, p1, 0.1)  
                add_farey_edge!(Graph_n, v_child3, p2, 0.1)
                update_edge_weight!(Graph_n, p1, p2, 1.0)  # Remove the old edge instead of setting weight to 0
            end 
        end
        return Graph_n, Value_n
    end
end

# ╔═╡ 244e477d-d8f6-487f-8b86-b59ecbebf3e2
FareyGraph(1)

# ╔═╡ 0768cb85-0f8b-4e42-9371-1ba1ea23069d
g = FareyGraph(2)[1]

# ╔═╡ 6a1423d2-887f-460f-8b3c-e9a8ff34ae6c
f = FareyGraph(3)[2]

# ╔═╡ 844a2b6b-b482-4f4b-9d2b-413d04ab0e4e
for e in edges(g)
	s = src(e)
	d = dst(e)
	println(first(f[s]))
end

# ╔═╡ a1a85fad-a318-4b81-aa18-1b071f4de67c
FareyGraph(2)

# ╔═╡ 4377aae9-ea7d-4644-89e1-b0c4536a0171
for e in edges(FareyGraph2(2)[1])
	println(e)
end


# ╔═╡ 081e2c92-5481-4359-b6db-b95462cbcfc0
@bind n PlutoUI.Slider(1:20)

# ╔═╡ e91d3115-e87f-4a3e-8f44-098512a9878b
FareyGraph2(2)[2]

# ╔═╡ ef42c906-c9d4-4266-b5be-265082cc8d1d
 Dict(k => string(v) for (k, v) in FareyGraph2(n)[2])

# ╔═╡ 04527d7f-9ef6-41e5-84eb-48f1b23a3382
gplot(FareyGraph2(n)[1])

# ╔═╡ ae1c785d-a677-4eb0-8d89-777fd8f3dba1
DictF = FareyGraph2(n)[2]

# ╔═╡ 3055b90b-f8fa-4b97-a9b9-3b9c0f5816b1
md"
## Farey Diagram angle functions : 
"

# ╔═╡ 271e24c9-62ba-4a19-8ceb-5c1343e6e00c
md"
### Normal angle implementation 
- The actual angles is a coefficient of the Farey equation 
"

# ╔═╡ f159435d-f940-4ffa-aafb-c443ec274ce7
function angle(a::Set{Rational{Int}})
	v = first(a)
	if v ∈ [Rational(-1,0),Rational(1,0)]
		return π
	elseif v>=-1 && v<=1
		return -float(v)*π/2
	else 
		return (float(1/v) + 2*sign(v)*1)*π/2 
		
	end
	
end

# ╔═╡ b414a07b-2f7a-4fe2-b83a-2ee34cd90ab2
Set(Rational(1,1))

# ╔═╡ d19011cd-4a4f-4502-8a56-9f88ee9f2cf7
md"
### Equidistant angle implementation
"

# ╔═╡ b941d450-82ff-47ab-a38b-d3e43ae059b7
md"""
Let $q \in \mathbb{Q}$ a fraction with continued fraction series $u_{n} = \{a_{0}:a_{1},a_{2},\dots, a_{n}\}$ :
Let $\Theta:\mathbb{Q} \Longrightarrow \mathbb{R}$ the angle of $q$ in the Farey Diagram :
$$\Theta(q) = \pi + \sum_{k= 1}^{n}\sum_{i = 1}^{a_{i}}(-1)^i(\frac{1}{2})^{k}\pi$$
( Note that the continued fraction series is finite because $q \in \mathbb{Q}$ )
"""

# ╔═╡ 9578f0da-d19e-4bc4-aba8-d70e680fc9ea
function continued_fraction(r::Rational)
    series = []
    num = numerator(r)
    den = denominator(r)
    
    while den != 0
        q = div(num, den)
        push!(series, q)
        num, den = den, num - q * den
    end
    
    return series
end

# ╔═╡ 6ef4e9ac-0867-4d44-9d62-3237417439e7
function theta(q) 
    cf = continued_fraction(q)
    angle = 0 # Start with π as per the formula
	if q == 0 
		return 0
	elseif q == 1 
		return π/2
	end
	a = []
	n = 2 
    for k in cf[2:end]
		a = cf[2:n]
		println(a)
		angle = angle +(-1)^(n)*((1/2)^sum(a))*π
		n = n+1
	end
	
	return angle
end

# ╔═╡ 0b7fde22-ed52-4648-8ee8-e02b3dd99066
function equi_angle(a::Set{Rational{Int}})
	v = first(a)
	if v ∈ [Rational(-1,0),Rational(1,0)]
		return π
	elseif v>=0 && v<=1
		return theta(v)
	elseif v>=-1 && v<0
		return -theta(-v)
	elseif v>1 
		return π - theta(1/v)
	else 
		return π + theta(-1/v)
	end
end

# ╔═╡ 2dfdbe63-02a0-471b-9147-fb0416d80d3e
equi_angle(Set([Rational(0,1)]))

# ╔═╡ 74e3475c-f443-417c-87e0-104a9cf02429
continued_fraction(1//1)

# ╔═╡ 789e2c22-59f0-45cf-9e7b-aa5ee0642e2b
π/2 - π/8

# ╔═╡ 1aa8a933-23bd-43ac-bec4-9b234383277b
theta(0//1)

# ╔═╡ c657bcac-2dd3-46ed-80a7-04ec00b420ea
π/2

# ╔═╡ e0220308-04d5-4936-b930-e6d400b2aa87


# ╔═╡ 9cbb0a91-66bc-4fb5-a5e3-539e148e0dde
for value in values(DictF)
	println(angle(value))
	
end 

# ╔═╡ b9efa191-cfa5-429b-8b1b-752286b96152
# ╠═╡ disabled = true
#=╠═╡

  ╠═╡ =#

# ╔═╡ 337898fd-a0a0-44fd-a311-fbc7d02e7cf6
my_radius = 400

# ╔═╡ 714d4bcd-1e57-4ac0-b0ab-69b3532d4e58
function tangent_line(x0, y0, h, k, r)
    A = x0 - h
    B = y0 - k
    C = r^2 + h*x0 + k*y0
    return A, B, -C
end

# ╔═╡ 313918c8-2867-483d-8945-56b4fb54243e
function line_intersection(A1, B1, C1, A2, B2, C2)
    M = [A1 B1; A2 B2]
    v = [C1; C2]
    return M \ v
end

# ╔═╡ 2a04392e-7edf-4b95-8a20-449db0d52a9f
function plot_farey_graph_with_edges2(m::Int,G)
    @svg begin 
    Drawing(1000, 1000)
    Luxor.origin()
    background("antiquewhite")
    DictG = G[2]
    graph = G[1]
    Edges= edges(graph)
    circle(Luxor.Point(0, 0), my_radius, action = :stroke)
    i = 1
    
    for e in Edges
        p1 = src(e)
        p2 = dst(e)
        x_1 = my_radius * cos(angle(DictG[p1]))
        y_1 = my_radius * sin(angle(DictG[p1]))
        x_2 = my_radius * cos(angle(DictG[p2]))
        y_2 = my_radius * sin(angle(DictG[p2]))
        angle1 = angle(DictG[p1])
        angle2 = angle(DictG[p2])
        coord1 = Luxor.Point(x_1,y_1)
        coord2 = Luxor.Point(x_2,y_2)
        if y_1== 0
            y_1 = 10^(-10)
        end
        ehehe_x = (y_1*(x_2^2) - y_2*(x_1^2)+(y_2-y_1)y_2*y_1)/(y_1*x_2 - y_2*x_1)

        ehehe_y = -(x_1/y_1)*ehehe_x +((x_1^2)/y_1) + y_1
		
        ehehe = Luxor.Point(ehehe_x,ehehe_y)
		println(ehehe)
        radius = distance(ehehe,coord1)
        f_1 = first(DictG[p1])
        f_2 = first(DictG[p2])
		#println(f_1,f_2,f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)],f_1 in [Rational(1, 0), Rational(-1, 0)],f_2 < 0,f_1 < 0)
		if (f_1,f_2) == (1,-1) || (f_1,f_2) ==(-1,1)
				Luxor.line(coord1,coord2, :stroke)
		
		elseif f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
			
                if f_1 in [Rational(1, 0), Rational(-1, 0)]
                    if f_2 < 0 
						
                        Luxor.arc2r(ehehe, coord1, coord2, :stroke)

                    else 

                        Luxor.arc2r(ehehe, coord2, coord1, :stroke)

                    end
                else
                    if f_1 < 0 

                        Luxor.arc2r(ehehe, coord2, coord1, :stroke)

                    else 

                        Luxor.arc2r(ehehe, coord1, coord2, :stroke)
                    end
                end

        elseif f_1< f_2
			
            Luxor.arc2r(ehehe, coord1, coord2, :stroke)
        else
			
            Luxor.arc2r(ehehe, coord2, coord1, :stroke)
        end
        #circle(ehehe,radius, :stroke)
		println("Drawing edge from $(f_1):$(coord1)  to $(f_2):$(coord2)")
		println("Center: $(ehehe),$(radius)")

    sethue("black")
    end
	for (k, pts) in DictG
        angle_rad = angle(pts)
        coord_x = my_radius * cos(angle_rad)
        coord_y = my_radius * sin(angle_rad)

        circle(Luxor.Point(coord_x, coord_y), 3, :fill)

        # Calculate text position
        text_radius = my_radius + 40
        coord_textx = text_radius * cos(angle_rad)
        coord_texty = text_radius * sin(angle_rad)

        # Create fraction text
        frac = first(pts)
        if denominator(frac) == 1
            label = string(numerator(frac))
        elseif frac ∈ [Rational(-1,0), Rational(1,0)]
            label = "∞"
        else
            label = string(numerator(frac), "/", denominator(frac))
        end

        # Add text
        fontsize(20)
        Luxor.text(label, Luxor.Point(coord_textx, coord_texty), halign=:center, valign=:middle)
         
    end
	finish()
end
end

# ╔═╡ 5c263bac-5fb6-49d6-b38d-2b067c4517e7
function plot_farey_graph_with_edges2equi(m::Int,G)
    @svg begin 
    Drawing(1000, 1000)
    Luxor.origin()
    background("antiquewhite")
    DictG = G[2]
    graph = G[1]
    Edges= edges(graph)
    circle(Luxor.Point(0, 0), my_radius, action = :stroke)
    i = 1
    
    for e in Edges
        p1 = src(e)
        p2 = dst(e)
        x_1 = my_radius * cos(equi_angle(DictG[p1]))
        y_1 = my_radius * sin(equi_angle(DictG[p1]))
        x_2 = my_radius * cos(equi_angle(DictG[p2]))
        y_2 = my_radius * sin(equi_angle(DictG[p2]))
        angle1 = equi_angle(DictG[p1])
        angle2 = equi_angle(DictG[p2])
        coord1 = Luxor.Point(x_1,y_1)
        coord2 = Luxor.Point(x_2,y_2)
        if y_1== 0
            y_1 = 10^(-10)
        end
        ehehe_x = (y_1*(x_2^2) - y_2*(x_1^2)+(y_2-y_1)y_2*y_1)/(y_1*x_2 - y_2*x_1)

        ehehe_y = -(x_1/y_1)*ehehe_x +((x_1^2)/y_1) + y_1
		
        ehehe = Luxor.Point(ehehe_x,ehehe_y)
		println(ehehe)
        radius = distance(ehehe,coord1)
        f_1 = first(DictG[p1])
        f_2 = first(DictG[p2])
		#println(f_1,f_2,f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)],f_1 in [Rational(1, 0), Rational(-1, 0)],f_2 < 0,f_1 < 0)
		if (f_1,f_2) == (1,-1) || (f_1,f_2) ==(-1,1)
				Luxor.line(coord1,coord2, :stroke)
		
		elseif f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
			
                if f_1 in [Rational(1, 0), Rational(-1, 0)]
                    if f_2 < 0 
						
                        Luxor.arc2r(ehehe, coord2, coord1, :stroke)

                    else 

                        Luxor.arc2r(ehehe, coord1, coord2, :stroke)

                    end
                else
                    if f_1 < 0 

                        Luxor.arc2r(ehehe, coord1, coord2, :stroke)

                    else 

                        Luxor.arc2r(ehehe, coord2, coord1, :stroke)
                    end
                end

        elseif f_1< f_2
			
            Luxor.arc2r(ehehe, coord2, coord1, :stroke)
        else
			
            Luxor.arc2r(ehehe, coord1, coord2, :stroke)
        end
        #circle(ehehe,radius, :stroke)
		println("Drawing edge from $(f_1):$(coord1)  to $(f_2):$(coord2)")
		println("Center: $(ehehe),$(radius)")

    sethue("black")
    end
	for (k, pts) in DictG
        angle_rad = equi_angle(pts)
        coord_x = my_radius * cos(angle_rad)
        coord_y = my_radius * sin(angle_rad)

        circle(Luxor.Point(coord_x, coord_y), 3, :fill)

        # Calculate text position
        text_radius = my_radius + 40
        coord_textx = text_radius * cos(angle_rad)
        coord_texty = text_radius * sin(angle_rad)

        # Create fraction text
        frac = first(pts)
        if denominator(frac) == 1
            label = string(numerator(frac))
        elseif frac ∈ [Rational(-1,0), Rational(1,0)]
            label = "∞"
        else
            label = string(numerator(frac), "/", denominator(frac))
        end

        # Add text
        fontsize(20)
        Luxor.text(label, Luxor.Point(coord_textx, coord_texty), halign=:center, valign=:middle)
         
    end
	finish()
end
end

# ╔═╡ 281c7482-bd09-488e-b28f-a086b8ba3810
function plot_farey_graph_with_edges3(m::Int,G)
    @svg begin 
    Drawing(1000, 1000)
    Luxor.origin()
    background("antiquewhite")
    DictG = G[2]
    graph = G[1]
    Edges= edges(graph)
    circle(Luxor.Point(0, 0), my_radius, action = :stroke)
    i = 1
    
    for e in Edges
        p1 = src(e)
        p2 = dst(e)
        x_1 = my_radius * cos(angle(DictG[p1]))
        y_1 = my_radius * sin(angle(DictG[p1]))
        x_2 = my_radius * cos(angle(DictG[p2]))
        y_2 = my_radius * sin(angle(DictG[p2]))
        angle1 = angle(DictG[p1])
        angle2 = angle(DictG[p2])
        coord1 = Luxor.Point(x_1,y_1)
        coord2 = Luxor.Point(x_2,y_2)
        if y_1== 0
            y_1 = 10^(-10)
        end
        ehehe_x = (y_1*(x_2^2) - y_2*(x_1^2)+(y_2-y_1)y_2*y_1)/(y_1*x_2 - y_2*x_1)

        ehehe_y = -(x_1/y_1)*ehehe_x +((x_1^2)/y_1) + y_1
		
        ehehe = Luxor.Point(ehehe_x,ehehe_y)
		println(ehehe)
        radius = distance(ehehe,coord1)
        f_1 = first(DictG[p1])
        f_2 = first(DictG[p2])
		#println(f_1,f_2,f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)],f_1 in [Rational(1, 0), Rational(-1, 0)],f_2 < 0,f_1 < 0)
		
		
        if f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
			
                if f_1 in [Rational(1, 0), Rational(-1, 0)]
                    if f_2 < 0 
						
                        Luxor.arc2r(ehehe, coord1, coord2, :stroke)

                    else 

                        Luxor.arc2r(ehehe, coord2, coord1, :stroke)

                    end
                else
                    if f_1 < 0 

                        Luxor.arc2r(ehehe, coord2, coord1, :stroke)

                    else 

                        Luxor.arc2r(ehehe, coord1, coord2, :stroke)
                    end
                end

        elseif f_1< f_2
			
            Luxor.arc2r(ehehe, coord1, coord2, :stroke)
        else
			
            Luxor.arc2r(ehehe, coord2, coord1, :stroke)
        end
        #circle(ehehe,radius, :stroke)
		println("Drawing edge from $(f_1) to $(f_2)")
		println("Coordinates: $(coord1) -> $(coord2)")
		println("Center point: $(ehehe)")
		println("Radius: $(radius)")
    sethue("black")
    end
	for (k, pts) in DictG
        angle_rad = angle(pts)
        coord_x = my_radius * cos(angle_rad)
        coord_y = my_radius * sin(angle_rad)

        circle(Luxor.Point(coord_x, coord_y), 3, :fill)

        # Calculate text position
        text_radius = my_radius + 40
        coord_textx = text_radius * cos(angle_rad)
        coord_texty = text_radius * sin(angle_rad)

        # Create fraction text
        frac = first(pts)
        if denominator(frac) == 1
            label = string(numerator(frac))
        elseif frac ∈ [Rational(-1,0), Rational(1,0)]
            label = "∞"
        else
            label = string(numerator(frac), "/", denominator(frac))
        end

        # Add text
        fontsize(20)
        Luxor.text(label, Luxor.Point(coord_textx, coord_texty), halign=:center, valign=:middle)
         
    end
	finish()
end
end

# ╔═╡ 16a0c42f-63da-40ae-9fc7-a2eccae6b34a
function plot_farey_graph_with_edges(m::Int)
	@svg begin 
    Drawing(1000, 1000)
    Luxor.origin()
    background("antiquewhite")
    DictG = FareyGraph(m)[2]
	graph = FareyGraph(m)[1]
    Edges= edges(graph)
    circle(Luxor.Point(0, 0), my_radius, action = :stroke)
    
    for (k, pts) in DictG
        angle_rad = angle(pts)
        coord_x = my_radius * cos(angle_rad)
        coord_y = my_radius * sin(angle_rad)
        
        circle(Luxor.Point(coord_x, coord_y), 3, :fill)
        
        # Calculate text position
        text_radius = my_radius + 40
        coord_textx = text_radius * cos(angle_rad)
        coord_texty = text_radius * sin(angle_rad)
        
        # Create fraction text
        frac = first(pts)
        if denominator(frac) == 1
            label = string(numerator(frac))
        elseif frac ∈ [Rational(-1,0), Rational(1,0)]
            label = "∞"
        else
            label = string(numerator(frac), "/", denominator(frac))
        end
        
        # Add text
        fontsize(20)
        Luxor.text(label, Luxor.Point(coord_textx, coord_texty), halign=:center, valign=:middle)
    end
	for e in Edges
		p1 = src(e)
		p2 = dst(e)
		x_1 = my_radius * cos(angle(DictG[p1]))
		y_1 = my_radius * sin(angle(DictG[p1]))
		x_2 = my_radius * cos(angle(DictG[p2]))
		y_2 = my_radius * sin(angle(DictG[p2]))
		coord1 = Luxor.Point(x_1,y_1)
        coord2 = Luxor.Point(x_2,y_2)
		ehehe_x = (y_1*(x_2^2) - y_2*(x_1^2))/(y_1*x_2 - y_2*x_1)
		ehehe_y =
		circle(coord1,coord2, :stroke)
		
	end
    
    finish()
end
end

# ╔═╡ 7575e497-e2c4-4526-8c88-f99509a7270a
@bind m PlutoUI.Slider(1:20)

# ╔═╡ 6456806b-34c5-4a10-a5d3-a0a7277164b9
plot_farey_graph_with_edges2(m,FareyGraph(m))

# ╔═╡ ea67c972-f01a-47b1-8e27-0ca81c71e839
plot_farey_graph_with_edges2equi(m,FareyGraph(m))

# ╔═╡ eb0c03ee-62c0-4cb9-b52f-1b1830eaa0a4
# ╠═╡ skip_as_script = true
#=╠═╡
function arc2r(scene, center, pt1, pt2; radius=nothing, linewidth=2, color=:black)
    if radius === nothing
        radius = norm(pt1 .- center)
    end
    
    # Calculate angles
    θ1 = atan(pt1[2] - center[2], pt1[1] - center[1])
    θ2 = atan(pt2[2] - center[2], pt2[1] - center[1])
    
    # Ensure the arc is drawn in the correct direction
    if θ2 < θ1
        θ2 += 2π
    end
    
    arc!(scene, Point2f0(center...), radius, θ1, θ2; linewidth, color)
end
  ╠═╡ =#

# ╔═╡ 292a52e9-7b90-4687-852d-b9eb5ad95fa3
# ╠═╡ skip_as_script = true
#=╠═╡
function plot_farey_graph(ax, m::Int, radius::Float64,n)
    
	if n == 1
    	DictG = FareyGraph(m)[2]
    	graph = FareyGraph(m)[1]
		c = :green
		
	else 
		DictG = FareyGraph2(m)[2]
		graph = FareyGraph2(m)[1]
		c = :blue
	end
    Edges = edges(graph)

    # Plot main circle
    θ = range(0, 2π, length=100)
    circle_x = radius .* cos.(θ)
    circle_y = radius .* sin.(θ)
	
	
	lines!(ax, circle_x, circle_y, color = :black)
	
    for (k, pts) in DictG
        angle_rad = angle(pts)
        coord_x = radius * cos(angle_rad)
        coord_y = radius * sin(angle_rad)

        GLMakie.scatter!(ax, [coord_x], [coord_y], color = :blue, markersize = 10)

        # Calculate text position
        text_radius = radius + 0.1 * radius
        coord_textx = text_radius * cos(angle_rad)
        coord_texty = text_radius * sin(angle_rad)

        # Create fraction text
        frac = first(pts)
        label = if denominator(frac) == 1
            string(numerator(frac))
        elseif frac ∈ [Rational(-1, 0), Rational(1, 0)]
            "∞"
        else
            string(numerator(frac), "/", denominator(frac))
        end

        # Add text
        text!(ax, label, position = (coord_textx, coord_texty), align = (:center, :center), fontsize = 14)
    end

    # Draw arcs for edges
    for e in Edges
        p1 = src(e)
        p2 = dst(e)
        x_1 = radius * cos(angle(DictG[p1]))
        y_1 = radius * sin(angle(DictG[p1]))
        x_2 = radius * cos(angle(DictG[p2]))
        y_2 = radius * sin(angle(DictG[p2]))
        coord1 = (x_1, y_1)
        coord2 = (x_2, y_2)

        if y_1 == 0
            y_1 = 1e-10
        end

        ehehe_x = (y_1 * (x_2^2) - y_2 * (x_1^2) + (y_2 - y_1) * y_2 * y_1) / (y_1 * x_2 - y_2 * x_1)
        ehehe_y = -(x_1 / y_1) * ehehe_x + ((x_1^2) / y_1) + y_1

        # Ensure valid center and radius calculations
        if !isfinite(ehehe_x) || !isfinite(ehehe_y)
            continue  # Skip this arc if the center is not valid
        end

        center = (ehehe_x, ehehe_y)
        circle_radius = sqrt((x_1 - ehehe_x)^2 + (y_1 - ehehe_y)^2)

        f_1 = first(DictG[p1])
        f_2 = first(DictG[p2])

        if f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
            if f_1 in [Rational(1, 0), Rational(-1, 0)]
                if f_2 < 0
                    arc2r(ax, center, coord1, coord2,color = c)
                else
                    arc2r(ax, center, coord2, coord1,color = c)
                end
            else
                if f_1 < 0
                    arc2r(ax, center, coord2, coord1,color = c)
                else
                    arc2r(ax, center, coord1, coord2,color = c)
                end
            end
        elseif f_1 < f_2
            arc2r(ax, center, coord1, coord2, color = c)
        else
            arc2r(ax, center, coord2, coord1, color = c)
        end

        # Draw the circle if valid
        
    end

    hidedecorations!(ax)
    hidespines!(ax)
    Makie.xlims!(ax, -radius * 1.2, radius * 1.2)
    Makie.ylims!(ax, -radius * 1.2, radius * 1.2)
end
  ╠═╡ =#

# ╔═╡ 3ac4dadf-28e1-4f08-8ebd-ebda389422b6
# ╠═╡ skip_as_script = true
#=╠═╡

function plot_farey_graph_equi(ax, m::Int, radius::Float64,n)
    
	if n == 1
    	DictG = FareyGraph(m)[2]
    	graph = FareyGraph(m)[1]
		c = :green
		
	else 
		DictG = FareyGraph2(m)[2]
		graph = FareyGraph2(m)[1]
		c = :blue
	end
    Edges = edges(graph)

    # Plot main circle
    θ = range(0, 2π, length=100)
    circle_x = radius .* cos.(θ)
    circle_y = radius .* sin.(θ)
	
	
	lines!(ax, circle_x, circle_y, color = :black)
	
    for (k, pts) in DictG
        angle_rad = equi_angle(pts)
        coord_x = radius * cos(angle_rad)
        coord_y = radius * sin(angle_rad)

        GLMakie.scatter!(ax, [coord_x], [coord_y], color = :blue, markersize = 10)

        # Calculate text position
        text_radius = radius + 0.1 * radius
        coord_textx = text_radius * cos(angle_rad)
        coord_texty = text_radius * sin(angle_rad)

        # Create fraction text
        frac = first(pts)
        label = if denominator(frac) == 1
            string(numerator(frac))
        elseif frac ∈ [Rational(-1, 0), Rational(1, 0)]
            "∞"
        else
            string(numerator(frac), "/", denominator(frac))
        end

        # Add text
        text!(ax, label, position = (coord_textx, coord_texty), align = (:center, :center), fontsize = 14)
    end

    # Draw arcs for edges
    for e in Edges
        p1 = src(e)
        p2 = dst(e)
        x_1 = radius * cos(equi_angle(DictG[p1]))
        y_1 = radius * sin(equi_angle(DictG[p1]))
        x_2 = radius * cos(equi_angle(DictG[p2]))
        y_2 = radius * sin(equi_angle(DictG[p2]))
        coord1 = (x_1, y_1)
        coord2 = (x_2, y_2)

        if y_1 == 0
            y_1 = 1e-10
        end

        ehehe_x = (y_1 * (x_2^2) - y_2 * (x_1^2) + (y_2 - y_1) * y_2 * y_1) / (y_1 * x_2 - y_2 * x_1)
        ehehe_y = -(x_1 / y_1) * ehehe_x + ((x_1^2) / y_1) + y_1

        # Ensure valid center and radius calculations
        if !isfinite(ehehe_x) || !isfinite(ehehe_y)
            continue  # Skip this arc if the center is not valid
        end

        center = (ehehe_x, ehehe_y)
        circle_radius = sqrt((x_1 - ehehe_x)^2 + (y_1 - ehehe_y)^2)

        f_1 = first(DictG[p1])
        f_2 = first(DictG[p2])

        if f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
            if f_1 in [Rational(1, 0), Rational(-1, 0)]
                if f_2 < 0
                    arc2r(ax, center, coord2, coord1,color = c)
                else
                    arc2r(ax, center, coord1, coord2,color = c)
                end
            else
                if f_1 < 0
                    arc2r(ax, center, coord1, coord2,color = c)
                else
                    arc2r(ax, center, coord2, coord1,color = c)
                end
            end
        elseif f_1 < f_2
            arc2r(ax, center, coord2, coord1, color = c)
        else
            arc2r(ax, center, coord1, coord2, color = c)
        end

        # Draw the circle if valid
        
    end

    hidedecorations!(ax)
    hidespines!(ax)
    Makie.xlims!(ax, -radius * 1.2, radius * 1.2)
    Makie.ylims!(ax, -radius * 1.2, radius * 1.2)
end
  ╠═╡ =#

# ╔═╡ 5f4f446f-306a-4094-b14d-86969a6fe111
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
function interactive_farey_graph()
    fig = Figure(size = (1920,1080))
    ax = Axis(fig[1, 1:3], aspect = DataAspect())  # Span the axis across more space
    
    m_slider = Makie.Slider(fig[2, 2], range = 1:20, startvalue = 1)
    n_slider = Makie.Slider(fig[3, 2], range = 1:2, startvalue = 1)
    both_toggle = Toggle(fig[2:3, 3], active = false)  # Toggle button on the side
    m_label = Label(fig[2, 1], "m:")
    n_label = Label(fig[3, 1], "n:")
    both_label = Label(fig[2:3, 4], "Show both n=1 and n=2:")
    
    function update_plot(m, n, both_active)
        empty!(ax)  # Clear the axis before plotting
        if both_active
            plot_farey_graph_equi(ax, m, 500.0, 1)
            plot_farey_graph_equi(ax, m, 500.0, 2)
        else
            plot_farey_graph_equi(ax, m, 500.0, n)
        end
    end
    
    on(m_slider.value) do m
        update_plot(m, n_slider.value[], both_toggle.active[])
    end
    
    on(n_slider.value) do n
        update_plot(m_slider.value[], n, both_toggle.active[])
    end
    
    on(both_toggle.active) do both_active
        update_plot(m_slider.value[], n_slider.value[], both_active)
    end
    
    # Initial plot
    update_plot(m_slider.value[], n_slider.value[], both_toggle.active[])
    
    # Enable zooming and panning
    display(fig)
end
  ╠═╡ =#

# ╔═╡ b1d69da2-0592-43f6-9734-06ee3d887311


# ╔═╡ d3bc6551-c75c-4ae7-822c-68d9cf45c73f


# ╔═╡ f329081b-6127-409c-811d-9a055d1ea1ce
# ╠═╡ disabled = true
#=╠═╡
m = 5 
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
GeometryTypes = "4d00f742-c7ba-57c2-abde-4428a4b178cb"
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
Luxor = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SimpleWeightedGraphs = "47aef6b3-ad0c-573a-a1e2-d07658019622"

[compat]
Colors = "~0.12.11"
GeometryTypes = "~0.8.5"
GraphPlot = "~0.5.2"
Graphs = "~1.9.0"
Luxor = "~4.1.0"
PlutoUI = "~0.7.59"
SimpleWeightedGraphs = "~1.4.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "8241f18bd755b4fa455bc3f672570221d6f624ee"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "bf6570a34c850f99407b494757f5d7ad233a7257"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.5"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.GeometryTypes]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "d796f7be0383b5416cd403420ce0af083b0f9b28"
uuid = "4d00f742-c7ba-57c2-abde-4428a4b178cb"
version = "0.8.5"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.GraphPlot]]
deps = ["ArnoldiMethod", "ColorTypes", "Colors", "Compose", "DelimitedFiles", "Graphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "5cd479730a0cb01f880eff119e9803c13f214cab"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.5.2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "899050ace26649433ef1af25bc17a815b3db52b7"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.9.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "ae0923dab7324e6bc980834f709c4cd83dd797ed"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.5+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Luxor]]
deps = ["Base64", "Cairo", "Colors", "DataStructures", "Dates", "FFMPEG", "FileIO", "PolygonAlgorithms", "PrecompileTools", "Random", "Rsvg"]
git-tree-sha1 = "134570038473304d709de27384621bd0810d23fa"
uuid = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
version = "4.1.0"

    [deps.Luxor.extensions]
    LuxorExtLatex = ["LaTeXStrings", "MathTeXEngine"]

    [deps.Luxor.weakdeps]
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    MathTeXEngine = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cb5a2ab6763464ae0f19c86c56c63d4a2b0f5bda"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.52.2+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PolygonAlgorithms]]
git-tree-sha1 = "a5ded6396172cff3bacdd1354d190b93cb667c4b"
uuid = "32a0d02f-32d9-4438-b5ed-3a2932b48f96"
version = "0.2.0"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rsvg]]
deps = ["Cairo", "Glib_jll", "Librsvg_jll"]
git-tree-sha1 = "3d3dc66eb46568fb3a5259034bfc752a0eb0c686"
uuid = "c4c386cf-5103-5370-be45-f3a111cca3b8"
version = "1.0.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "4b33e0e081a825dbfaf314decf58fa47e53d6acb"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.4.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "d9717ce3518dc68a99e6b96300813760d887a01d"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.1+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "86e7731be08b12fa5e741f719603ae740e16b666"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.10+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╠═15f0e8a0-59e4-11ef-288b-17fd5d6af339
# ╠═42a780fe-1659-4d9b-b900-3321ce025dd6
# ╠═0dfb78a9-a549-4c56-92f5-f4e5aa6f26d5
# ╠═68057bda-5823-4e28-8a56-cfdfedb15e5d
# ╠═70d7cb70-f755-4241-bbbd-77852f0b99d5
# ╠═817a2afb-3845-4cc7-bf9d-cb6f78b50de6
# ╠═ded767d6-493f-40af-bba6-7f3bd5f9b57a
# ╠═2c41b5bc-794a-432e-8076-5be96de4bc70
# ╠═b0985351-2382-4ecb-9bb1-1ae8bb9271a9
# ╟─a6989575-1ced-42ec-bb25-dcde4a750b71
# ╠═8183a1a1-af93-4e48-bf57-7d1443a9e79e
# ╠═d5de3fdc-ec9d-40ac-8ba2-ff1db06e8c72
# ╠═244e477d-d8f6-487f-8b86-b59ecbebf3e2
# ╠═0768cb85-0f8b-4e42-9371-1ba1ea23069d
# ╠═6a1423d2-887f-460f-8b3c-e9a8ff34ae6c
# ╠═844a2b6b-b482-4f4b-9d2b-413d04ab0e4e
# ╠═a1a85fad-a318-4b81-aa18-1b071f4de67c
# ╠═4377aae9-ea7d-4644-89e1-b0c4536a0171
# ╠═081e2c92-5481-4359-b6db-b95462cbcfc0
# ╠═e91d3115-e87f-4a3e-8f44-098512a9878b
# ╠═ef42c906-c9d4-4266-b5be-265082cc8d1d
# ╠═04527d7f-9ef6-41e5-84eb-48f1b23a3382
# ╠═ae1c785d-a677-4eb0-8d89-777fd8f3dba1
# ╠═3055b90b-f8fa-4b97-a9b9-3b9c0f5816b1
# ╠═271e24c9-62ba-4a19-8ceb-5c1343e6e00c
# ╠═0b7fde22-ed52-4648-8ee8-e02b3dd99066
# ╠═2dfdbe63-02a0-471b-9147-fb0416d80d3e
# ╠═f159435d-f940-4ffa-aafb-c443ec274ce7
# ╠═b414a07b-2f7a-4fe2-b83a-2ee34cd90ab2
# ╠═d19011cd-4a4f-4502-8a56-9f88ee9f2cf7
# ╠═b941d450-82ff-47ab-a38b-d3e43ae059b7
# ╠═6ef4e9ac-0867-4d44-9d62-3237417439e7
# ╠═9578f0da-d19e-4bc4-aba8-d70e680fc9ea
# ╠═74e3475c-f443-417c-87e0-104a9cf02429
# ╠═789e2c22-59f0-45cf-9e7b-aa5ee0642e2b
# ╠═1aa8a933-23bd-43ac-bec4-9b234383277b
# ╠═c657bcac-2dd3-46ed-80a7-04ec00b420ea
# ╠═e0220308-04d5-4936-b930-e6d400b2aa87
# ╠═9cbb0a91-66bc-4fb5-a5e3-539e148e0dde
# ╠═b9efa191-cfa5-429b-8b1b-752286b96152
# ╠═337898fd-a0a0-44fd-a311-fbc7d02e7cf6
# ╠═714d4bcd-1e57-4ac0-b0ab-69b3532d4e58
# ╠═313918c8-2867-483d-8945-56b4fb54243e
# ╠═2a04392e-7edf-4b95-8a20-449db0d52a9f
# ╠═5c263bac-5fb6-49d6-b38d-2b067c4517e7
# ╟─281c7482-bd09-488e-b28f-a086b8ba3810
# ╟─16a0c42f-63da-40ae-9fc7-a2eccae6b34a
# ╠═7575e497-e2c4-4526-8c88-f99509a7270a
# ╠═6456806b-34c5-4a10-a5d3-a0a7277164b9
# ╠═ea67c972-f01a-47b1-8e27-0ca81c71e839
# ╟─b03a4559-a960-4b56-97ed-a144407ba2eb
# ╠═eb0c03ee-62c0-4cb9-b52f-1b1830eaa0a4
# ╠═292a52e9-7b90-4687-852d-b9eb5ad95fa3
# ╠═3ac4dadf-28e1-4f08-8ebd-ebda389422b6
# ╠═5f4f446f-306a-4094-b14d-86969a6fe111
# ╠═b1d69da2-0592-43f6-9734-06ee3d887311
# ╠═d3bc6551-c75c-4ae7-822c-68d9cf45c73f
# ╠═f329081b-6127-409c-811d-9a055d1ea1ce
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
