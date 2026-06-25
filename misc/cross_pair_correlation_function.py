import numpy as np
from scipy.spatial.distance import cdist
from muspan.query.get_centroids import get_centroids
from muspan.helpers.object_to_object_distance_matrix import object_to_object_distance_matrix
from muspan.spatial_statistics.helpers.update_circle_areas_using_boundaries import update_circle_areas_using_boundaries
from muspan.query.return_object_IDs_from_query_like import return_object_IDs_from_query_like
from muspan.query.filter_objects_to_included_regions import filter_objects_to_included_regions
from muspan.helpers.combine_include_and_exclude_boundaries import combine_include_and_exclude_boundaries
from muspan.helpers.clean_up import clean_up
from muspan.visualise.visualise_spatial_statistic import visualise_spatial_statistic
from muspan.region_based.generate_quadrats import generate_quadrats
from muspan.spatial_statistics.helpers import get_shape_boundaries

def cross_pair_correlation_function(domain,population_A,population_B,include_boundaries=None,exclude_boundaries=None,boundary_exclude_distance=0,distance_metric='euclidean',max_R=100,annulus_step=10,annulus_width=10,exclude_zero=False,remain_within_connected_component=False, return_confidence_interval=False,return_PCF_contributions=False,visualise_output=False,visualise_spatial_statistic_kwargs=dict()):
    
    object_indices_A_unfiltered = return_object_IDs_from_query_like(domain,population_A)
    object_indices_B_unfiltered = return_object_IDs_from_query_like(domain,population_B)
    include_region_indices = combine_include_and_exclude_boundaries(domain,include_boundaries,exclude_boundaries)

    object_indices_A = filter_objects_to_included_regions(domain, object_indices_A_unfiltered, include_region_indices, boundary_exclude_distance=boundary_exclude_distance)
    object_indices_B = filter_objects_to_included_regions(domain, object_indices_B_unfiltered, include_region_indices, boundary_exclude_distance=boundary_exclude_distance)
    
    p_A, inds_A = get_centroids(domain, object_indices_A)
    p_B, inds_B = get_centroids(domain, object_indices_B)
    
    total_include_analysis_area = 0
    for i in include_region_indices:
        total_include_analysis_area+=domain.objects[i].area
    # Get distance of each point to the include_region boundaries
    distances_to_boundaries = object_to_object_distance_matrix(domain,inds_A,include_region_indices)

    # Sort out the inner and outer radii of the annulus
    radii = np.arange(0, max_R+annulus_step, annulus_step)
    radii_inner = radii - annulus_width/2
    radii_outer = radii + annulus_width/2

    # Get areas of circles around each point within these boundaries, accounting for holes
    # Unique radii - any radius considered in either inner or outer
    unique_radii = np.unique(np.hstack((radii_inner,radii_outer)))
    #if 0 in unique_radii:
    #    unique_radii = unique_radii[1:]
    unique_radii = unique_radii[unique_radii > 0]
    
    areas = np.ones(shape=(len(inds_A),len(unique_radii)))
    
    # Default area is for a circle of radius r_i
    areas = areas*np.pi*unique_radii**2
    
    # Now boundary correct relevant areas
    boundaries_as_shapes = get_shape_boundaries(domain,include_region_indices)
    for r_i, r in enumerate(unique_radii):   
        areas_temp = areas[:,r_i]
        areas_A=update_circle_areas_using_boundaries(p_A,areas_temp,distances_to_boundaries,r,remain_within_connected_component,boundaries_as_shapes)
        areas[:,r_i] = areas_A
     
        
    # Now do PCF
    density_B = np.shape(p_B)[0] / total_include_analysis_area
    if remain_within_connected_component:
        #TODO implement this - check containing component of all points, only include pairs which are within the same one
        raise NotImplementedError('remain_within_connected_component is currently not implemented')
    distances_A_to_B = cdist(p_A, p_B, metric=distance_metric)
    # Don't measure from an object to itself
    for a_ind, a in enumerate(inds_A):
        if a in inds_B:
            b_ind = np.argwhere(inds_B == a)[0][0]
            # Set this to be well outside the final annulus
            distances_A_to_B[a_ind,b_ind] = radii_outer[-1]+max_R
            
            

    cross_PCF_A_to_B = np.zeros(shape=(len(radii_outer),1))

    contributions = np.zeros(shape=(len(inds_A),len(radii_outer)))
    for annulus in range(len(radii_outer)):
        # print(annulus)
        n_pts = len(inds_A)
        inner = radii_inner[annulus]
        outer = radii_outer[annulus]
        distanceMask = np.logical_and((distances_A_to_B >= inner),(distances_A_to_B < outer))
        if outer <= 0 and exclude_zero:
            distanceMask[distances_A_to_B == 0] = False
        # Get indices of which unique_radii the inner and outer annuli correspond to
        if inner > 0 and outer >0:
            #print(inner)
            #print(unique_radii)
            inner_index = np.where(unique_radii == inner)[0][0]
            outer_index = np.where(unique_radii == outer)[0][0]
            assert(outer_index > inner_index)
            # Get areas of annulii
            annulus_areas = areas[:,outer_index] - areas[:,inner_index]
        elif inner <= 0 and outer > 0:
            
            outer_index = np.where(unique_radii == outer)[0][0]
            # Get areas of annulii
            annulus_areas = areas[:,outer_index]
        else:
            annulus_areas = np.zeros(len(inds_A))
     
        
        zero_indices = np.where(annulus_areas == 0)[0]
        # It's possible that some annulus_areas are 0 - i.e., for that point, none of the domain is within the annulus
        # Replace those values with 1 so that we don't get an error here...
        if len(zero_indices) > 0:
            n_pts = n_pts - len(zero_indices)
            if n_pts == 0:
                print(f'A PCF annulus is entirely outside the calculation area (outer radius = {radii_outer[annulus]}); consider reducing max_R')
            for z in zero_indices:
                annulus_areas[z] = 1

        contribution = np.sum(distanceMask, axis=1) / (density_B*annulus_areas)
        # ...ensure those values don't add to the contribution to the cross PCF...
        if len(zero_indices)>0:
            for z in zero_indices:
                contribution[z] = 0
        cross_PCF_A_to_B[annulus] += np.sum(contribution, axis=0)/n_pts
        if return_PCF_contributions or return_confidence_interval:
            # ...and finally set them to nan if we're returning contributions
            if len(zero_indices)>0:
                contribution[zero_indices] = 0
            contributions[:,annulus] = contribution
    # cross_PCF_A_to_B /= len(inds_A)
    
    
    
    if return_confidence_interval:
        # Block bootstrapping - see Loh et al 2008
        #TODO in future this could look at a query to decide which points are in which block
        # For now, use square lattice with edge length 20% of the domain side (5x5 lattice)
        region_kwargs = {}
        region_kwargs['include_boundaries']=include_boundaries
        region_kwargs['exclude_boundaries']=exclude_boundaries
        region_kwargs['population']=object_indices_A
        region_kwargs['regions_collection_name']='PCF bootstrap regions collection'
        region_kwargs['return_IDs']=True
        region_kwargs['remove_empty_regions']=False
        region_kwargs['region_label_name']='PCF bootstrap regions'
        region_kwargs['region_include_method']='clip'
        domain_short_side = np.min((domain.bounding_box[1,0]-domain.bounding_box[0,0],domain.bounding_box[1,1]-domain.bounding_box[0,1]))
        region_kwargs['side_length']=0.2*domain_short_side
        regions_object_ids=generate_quadrats(domain,assign_objects_using_labels=True,**region_kwargs)
        
        from muspan.query.get_labels import get_labels
        labs, OIs = get_labels(domain, 'PCF bootstrap regions')
        mask = [v in object_indices_A for v in OIs]
        OIs = OIs[mask]
        labs = labs[mask]
        
        region_options = domain.labels['PCF bootstrap regions']['categories']
        nQuadrats = len(region_options) # How many quadrats do we need to sample per PCF?
        numBootstrapSims = 1000
        samplePCFs = np.zeros(shape=(numBootstrapSims, np.shape(contributions)[1]))
        toSampleStr = np.random.choice(region_options,size=(nQuadrats,numBootstrapSims))
        toSample = np.zeros(shape=(nQuadrats,numBootstrapSims),dtype=int)
        
        quad_contributions = np.zeros((nQuadrats,np.shape(contributions)[1]))
        quadNs = np.zeros(nQuadrats)
        for j in range(len(region_options)):
            quadID = region_options[j]
            sampleMask = toSampleStr == quadID
            toSample[sampleMask] = j
            
            accept = labs==quadID
            if sum(accept) > 0:
                quad_contributions[j,:] = np.sum(contributions[accept,:],axis=0)
                quadNs[j] = sum(accept)
        sample = np.sum(quad_contributions[toSample,:],axis=0)
        Ns = np.sum(quadNs[toSample],axis=0)

        samplePCFs = sample / Ns[:,np.newaxis]

        # Get 95% CI
        PCF_min = 2*cross_PCF_A_to_B.T - np.percentile(samplePCFs, 97.5, axis=0)
        PCF_max = 2*cross_PCF_A_to_B.T - np.percentile(samplePCFs, 2.5, axis=0)
        PCF_min[PCF_min<0] = 0
        confidence_intervals = np.array([PCF_min[0], PCF_max[0]])
        
        # Clean up the domain
        domain.delete_objects(('Collection','PCF bootstrap regions collection'))
        domain.delete_labels('PCF bootstrap regions')    
    
    # Finally, delete the collection '__temporary_boundary_shapes'
    clean_up(domain)
    

    
    # visualise the results
    if visualise_output:
        if not return_confidence_interval:
            confidence_intervals = None
        visualise_spatial_statistic_kwargs=dict(confidence_intervals=confidence_intervals)|visualise_spatial_statistic_kwargs
        visualise_spatial_statistic(radii,cross_PCF_A_to_B,statistic_to_plot='PCF',**visualise_spatial_statistic_kwargs)
    # There must be a better way of doing this
    if (not return_PCF_contributions) and (not return_confidence_interval):
        # Standard return
        return radii, cross_PCF_A_to_B.flatten()
    else:
        # We're returning at least one of the two options
        if return_PCF_contributions and (not return_confidence_interval):
            return radii, cross_PCF_A_to_B.flatten(), contributions
        elif (not return_PCF_contributions) and return_confidence_interval:
            return radii, cross_PCF_A_to_B.flatten(), confidence_intervals
        else:
            # Return everything
            assert(return_PCF_contributions and return_confidence_interval)
            return radii, cross_PCF_A_to_B.flatten(), confidence_intervals, contributions