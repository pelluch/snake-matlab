function x  = CrossDetector(P)
% simple_Polygon(): test if a Polygon is simple or not
%     Input:  Pn = a polygon with n vertices V[]
%     Return: FALSE(0) = is NOT simple
%             TRUE(1)  = IS simple
EventQueue  Eq(Pn);
SweepLine   SL(Pn);
Event*      e; %  the current event
SLseg*      s; %  the current SL segment

% This loop processes all events in the sorted queue
% Events are only left or right vertices since
% No new events will be added (an intersect => Done)

while (e = Eq.next())          % while there are events
    if (e.type == LEFT)       % process a left vertex
        s = SL.add(e);          % add it to the sweep line
        if (SL.intersect(  s, s.above))
            return FALSE;      % Pn is NOT simple
        end
        if (SL.intersect(  s, s.below))
            return FALSE;      % Pn is NOT simple
        end
        
        
    else                       % processs a right vertex
        s = SL.find(e);
        if (SL.intersect(  s.above, s.below))
            return FALSE;      % Pn is NOT simple
            SL.remove(s);           % remove it from the sweep line
        end
    end
end


return TRUE;      % Pn IS simple




% EventQueue Class

% Event element data struct
typedef struct _event Event;
struct _event {
    int      edge;          % polygon edge i is V[i] to V[i+1]
    int      type;          % event type: LEFT or RIGHT vertex
    Point*   eV;            % event vertex
};

int E_compare( const void* v1, const void* v2 ) % qsort compare two events
{
    Event**    pe1 = (Event**)v1;
    Event**    pe2 = (Event**)v2;

    return xyorder( (*pe1).eV, (*pe2).eV );
}

% the EventQueue is a presorted array (no insertions needed)
class EventQueue {
    int      ne;                % total number of events in array
    int      ix;                % index of next event on queue
    Event*   Edata;             % array of all events
    Event**  Eq;                % sorted list of event pointers
public:
              EventQueue(Polygon P);     % constructor
             ~EventQueue(void)           % destructor
                  { delete[] Eq; delete[] Edata;}

    Event*   next();                     % next event on queue
};

% EventQueue Routines
EventQueue::EventQueue( Polygon P )
{
    ix = 0;
    ne = 2 * P.n;           % 2 vertex events for each edge
    Edata = (Event*)new Event[ne];
    Eq = (Event**)new (Event*)[ne];
    for (int i=0; i < ne; i++)           % init Eq array pointers
        Eq[i] = &Edata[i];

    % Initialize event queue with edge segment endpoints
    for (int i=0; i < P.n; i++) {        % init data for edge i
        Eq[2*i].edge = i;
        Eq[2*i+1].edge = i;
        Eq[2*i].eV   = &(P.V[i]);
        Eq[2*i+1].eV = &(P.V[i+1]);
        if (xyorder( &P.V[i], &P.V[i+1]) < 0)  { % determine type
            Eq[2*i].type    = LEFT;
             Eq[2*i+1].type = RIGHT;
        }
        else {
            Eq[2*i].type    = RIGHT;
             Eq[2*i+1].type = LEFT;
        }
    }
    % Sort Eq[] by increasing x and y
    qsort( Eq, ne, sizeof(Event*), E_compare );
}

Event* EventQueue::next()
{
    if (ix >= ne)
        return (Event*)0;
    else
        return Eq[ix++];
}
%===================================================================
 


% SweepLine Class

% SweepLine segment data struct
typedef struct _SL_segment SLseg;
struct _SL_segment {
    int      edge;          % polygon edge i is V[i] to V[i+1]
    Point    lP;            % leftmost vertex point
    Point    rP;            % rightmost vertex point
    SLseg*   above;         % segment above this one
    SLseg*   below;         % segment below this one
};

% the Sweep Line itself
class SweepLine {
    int      nv;            % number of vertices in polygon
    Polygon* Pn;            % initial Polygon
    BBT      Tree;          % balanced binary tree
public:
              SweepLine(Polygon P)            % constructor
                  { nv = P.n; Pn = &P; }
             ~SweepLine(void)                 % destructor
                  { Tree.freetree();}

    SLseg*   add( Event* );
    SLseg*   find( Event* );
    int      intersect( SLseg*, SLseg*  );
    void     remove( SLseg* );
};

SLseg* SweepLine::add( Event* E )
{
    % fill in SLseg element data
    SLseg* s = new SLseg;
    s.edge  = E.edge;

    % if it is being added, then it must be a LEFT edge event
    % but need to determine which endpoint is the left one 
    Point* v1 = &(Pn.V[s.edge]); 
    Point* v2 = &(Pn.V[s.edge+1]); 
    if (xyorder( v1, v2) < 0) { % determine which is leftmost
        s.lP = *v1; 
        s.rP = *v2; 
    } 
    else { 
        s.rP = *v1; 
        s.lP = *v2;
    }
    s.above = (SLseg*)0;
    s.below = (SLseg*)0;

    % add a node to the balanced binary tree
    Tnode* nd = Tree.insert(s);
    Tnode* nx = Tree.next(nd);
    Tnode* np = Tree.prev(nd);
    if (nx != (Tnode*)0) {
        s.above = (SLseg*)nx.val;
        s.above.below = s;
    }
    if (np != (Tnode*)0) {
        s.below = (SLseg*)np.val;
        s.below.above = s;
    }
    return s;
}

function y = find( E )

% need a segment to find it in the tree
s.edge  = E.edge;
s.above = 0;
s.below = 0;

nd = find(s);

if (nd == 0)
    y = 0;
    return;
end

y = nd.val;

return;

function y = remove(  s )
{
    % remove the node from the balanced binary tree
    nd = find(s);
    if (nd == 0)
        return;       % not there

    % get the above and below segments pointing to each other
    Tnode* nx = Tree.next(nd);
    if (nx != (Tnode*)0) {
        SLseg* sx = (SLseg*)(nx.val);
        sx.below = s.below;
    }
    Tnode* np = Tree.prev(nd);
    if (np != (Tnode*)0) {
        SLseg* sp = (SLseg*)(np.val);
        sp.above = s.above;
    }
    Tree.remove(nd);       % now  can safely remove it
    delete s;
}

% test intersect of 2 segments and return: 0=none, 1=intersect


function y = Intersect( s1, s2)

    if (s1 == [] || s2 == [])
        y = false;
        return;       % no intersect if either segment doesn't exist
    end

    % check for consecutive edges in polygon
e1 = s1.edge;
e2 = s2.edge;
if ((mod(e1+1,nv) == e2) || (e1 == mod(e2+1,nv)))
    y = false;
    return;       % no non-simple intersect since consecutive
end

    % test for existence of an intersect point
    
    lsign = crosspod(s1.lP, s1.rP, s2.lP);    %  s2 left point sign
    rsign = crosspod(s1.lP, s1.rP, s2.rP);    %  s2 right point sign
    if (lsign * rsign > 0) % s2 endpoints have same sign  relative to s1
        y = false;
        return;       % => on same side => no intersect is possible
    end
    
    lsign = crosspod(s2.lP, s2.rP, s1.lP);    %  s1 left point sign
    rsign = crosspod(s2.lP, s2.rP, s1.rP);    %  s1 right point sign
    if (lsign * rsign > 0) % s1 endpoints have same sign  relative to s2
        y = false;
        return;       % => on same side => no intersect is possible
    end
    % the segments s1 and s2 straddle each other
    y = true;         % => an intersect exists
end
    
%===================================================================
 

