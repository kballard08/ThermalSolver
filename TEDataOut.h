/*
 * TMDataOut.h
 *
 *  Created on: Apr 23, 2013
 *      Author: kballard
 */

#ifndef TEDATAOUT_H_
#define TEDATAOUT_H_

#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/data_component_interpretation.h>

#include <deal.II/base/std_cxx1x/shared_ptr.h>

#include "TEPostProcessor.h"

namespace FEASolverNS
{

using namespace dealii;

/*
 * The point of this class is to allow the compute_derived_quanitites in the
 * Postprocessor to know the material_id of the cell on which it is operating.
 * To do this, the TEDataOut class derives from the DataOut class and overrides
 * the build patch methods. Only build_one_patch needed to be overriden, but
 * by the design of DataOut the method and a few other needed methods were private
 * in scope.  So The code was copied to this class.  I (Keith Ballard) am not the
 * author of the append_patch_to_list function or the build_patches and iterator
 * methods of this class.  I copied them from DataOut because I didn't have a
 * choice given DataOut's design.  I should point out that I heavily modified
 * build_one_patch but it too is based off teh DataOut private method.
 *
 * The TEDataOut class ONLY works with the TEPostProcessor class as its
 * PostProcessor.  Any other class or the absense of it will throw an exception
 * with a message explaining what to change.  This is hackish I admit, but I
 * can't find a better solution without changing the actual source code of deal ii,
 * which I really don't want to do (due to compatibility across machines where
 * it is already installed).
 *
 * This class is not very generic, nor is it meant to.  It is valid with a specific
 * Postprocessor and requires vector outputs (scalar will throw an error).  Hopefully,
 * the DataOut class can be modified in the future to accomodate polymorphism on its
 * currently private methods.  Also, append_patch_to_list really needs to move to the
 * header file rather than the source file in deal II, that is why I had to copy it
 * here.
 */

// Declare append_patch_to_list
// Definition later
template <int dim, int spacedim>
void append_patch_to_list (const DataOutBase::Patch<dim,spacedim> &patch,
					  std::vector<DataOutBase::Patch<dim,spacedim> > &patches);

template <int dim, class DH=DoFHandler<dim> >
class TEDataOut : public DataOut<dim, DH>
{
public:
	typedef typename DataOut<dim, DH>::cell_iterator cell_iterator;
	typedef typename DataOut<dim, DH>::active_cell_iterator active_cell_iterator;
	typedef typename DataOut<dim, DH>::CurvedCellRegion CurvedCellRegion;

	virtual void build_patches (const unsigned int n_subdivisions = 0);
	virtual void build_patches (const Mapping<DH::dimension,DH::space_dimension> &mapping,
	                              const unsigned int n_subdivisions = 0,
	                              const CurvedCellRegion curved_region = DataOut<dim, DH>::curved_boundary);

	/**
	* Exception
	*/
	DeclException1 (ExcInvalidNumberOfSubdivisions,
				  int,
				  << "The number of subdivisions per patch, " << arg1
				  << ", is not valid.");
private:
	/**
	* Return the first cell produced
	* by the
	* first_cell()/next_cell()
	* function pair that is locally
	* owned. If this object operates
	* on a non-distributed
	* triangulation, the result
	* equals what first_cell()
	* returns.
	*/
	cell_iterator first_locally_owned_cell ();

	/**
	* Return the next cell produced
	* by the next_cell() function
	* that is locally owned. If this
	* object operates on a
	* non-distributed triangulation,
	* the result equals what
	* first_cell() returns.
	*/
	cell_iterator next_locally_owned_cell (const cell_iterator &cell);

	void build_one_patch (const std::pair<cell_iterator, unsigned int> *cell_and_index,
						internal::DataOut::ParallelData<DH::dimension, DH::space_dimension> &data,
						::dealii::DataOutBase::Patch<DH::dimension, DH::space_dimension> &patch,
						const CurvedCellRegion curved_cell_region);
};

template <int dim, class DH>
void TEDataOut<dim,DH>::build_patches (const unsigned int n_subdivisions) {
	build_patches (StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping,
	                 n_subdivisions, DataOut<dim, DH>::no_curved_cells);
}

template <int dim, class DH>
void TEDataOut<dim,DH>::build_patches (const Mapping<DH::dimension,DH::space_dimension> &mapping,
                                     const unsigned int nnnn_subdivisions,
                                     const CurvedCellRegion curved_region)
{
	// Check consistency of redundant
	// template parameter
	Assert (dim==DH::dimension, ExcDimensionMismatch(dim, DH::dimension));

	typedef DataOut_DoFData<DH, DH::dimension, DH::space_dimension> BaseClass;
	Assert (this->dofs != 0, typename BaseClass::ExcNoDoFHandlerSelected());

	const unsigned int n_subdivisions = (nnnn_subdivisions != 0)
									  ? nnnn_subdivisions
									  : this->default_subdivisions;
	Assert (n_subdivisions >= 1,
		  ExcInvalidNumberOfSubdivisions(n_subdivisions));

	// first count the cells we want to
	// create patches of. also fill the
	// object that maps the cell
	// indices to the patch numbers, as
	// this will be needed for
	// generation of neighborship
	// information
	std::vector<std::vector<unsigned int> > cell_to_patch_index_map;
	cell_to_patch_index_map.resize (this->dofs->get_tria().n_levels());
	for (unsigned int l=0; l<this->dofs->get_tria().n_levels(); ++l)
    {
		unsigned int max_index = 0;
		for (cell_iterator cell=first_locally_owned_cell(); cell != this->dofs->end();
		   cell = next_locally_owned_cell(cell))
		if (static_cast<unsigned int>(cell->level()) == l)
		  max_index = std::max (max_index,
								static_cast<unsigned int>(cell->index()));

		cell_to_patch_index_map[l].resize (max_index+1,
										 dealii::DataOutBase::Patch<DH::dimension,DH::space_dimension>::no_neighbor);
    }

	std::vector<std::pair<cell_iterator, unsigned int> > all_cells;
	{
		// set the index of the first
		// cell. if
		// first_locally_owned_cell /
		// next_locally_owned_cell
		// returns non-active cells, then
		// the index is not usable
		// anyway, but otherwise we
		// should keep track where we are
		unsigned int index;
		if ((first_locally_owned_cell() == this->dofs->end())
			||
			(first_locally_owned_cell()->has_children()))
		  index = 0;
		else
		  index = std::distance (this->dofs->begin_active(),
								 active_cell_iterator(first_locally_owned_cell()));
		for (cell_iterator cell=first_locally_owned_cell(); cell != this->dofs->end();
			 cell = next_locally_owned_cell(cell))
		{
			Assert (static_cast<unsigned int>(cell->level()) <
					cell_to_patch_index_map.size(),
					ExcInternalError());
			Assert (static_cast<unsigned int>(cell->index()) <
					cell_to_patch_index_map[cell->level()].size(),
					ExcInternalError());

			cell_to_patch_index_map[cell->level()][cell->index()] = all_cells.size();

			all_cells.push_back (std::make_pair(cell, index));

			// if both this and the next
			// cell are active, then
			// increment the index that
			// keeps track on which
			// active cell we are sitting
			// correctly. if one of the
			// cells is not active, then
			// this index doesn't mean
			// anything anyway, so just
			// ignore it. same if we are
			// at the end of the range
			if (!cell->has_children() &&
				next_locally_owned_cell(cell) != this->dofs->end() &&
				!next_locally_owned_cell(cell)->has_children())
			  index += std::distance (active_cell_iterator(cell),
								  	  active_cell_iterator(next_locally_owned_cell(cell)));
		}
	}

	this->patches.clear ();
	this->patches.reserve (all_cells.size());
	Assert (this->patches.size() == 0, ExcInternalError());

	// now create a default patch and a
	// default object for the
	// WorkStream object to work with
	const QTrapez<1>     q_trapez;
	const QIterated<DH::dimension> patch_points (q_trapez, n_subdivisions);

	const unsigned int n_components   = this->dofs->get_fe().n_components();
	unsigned int n_datasets=this->cell_data.size();
	for (unsigned int i=0; i<this->dof_data.size(); ++i)
	n_datasets += this->dof_data[i]->n_output_variables;

	std::vector<unsigned int> n_postprocessor_outputs (this->dof_data.size());
	for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
	if (this->dof_data[dataset]->postprocessor)
	  n_postprocessor_outputs[dataset] = this->dof_data[dataset]->n_output_variables;
	else
	  n_postprocessor_outputs[dataset] = 0;

	const CurvedCellRegion curved_cell_region
	= (n_subdivisions<2 ? DataOut<dim, DH>::no_curved_cells : curved_region);

	UpdateFlags update_flags = update_values;
	if (curved_cell_region != DataOut<dim, DH>::no_curved_cells)
	update_flags |= update_quadrature_points;

	for (unsigned int i=0; i<this->dof_data.size(); ++i)
	if (this->dof_data[i]->postprocessor)
	  update_flags |= this->dof_data[i]->postprocessor->get_needed_update_flags();
	// perhaps update_normal_vectors is present,
	// which would only be useful on faces, but
	// we may not use it here.
	Assert (!(update_flags & update_normal_vectors),
		  ExcMessage("The update of normal vectors may not be requested for evaluation of "
					 "data on cells via DataPostprocessor."));

	internal::DataOut::ParallelData<DH::dimension, DH::space_dimension>
	thread_data (patch_points,
			   n_components, n_datasets, n_subdivisions,
			   n_postprocessor_outputs,
			   mapping,
			   cell_to_patch_index_map,
			   this->dofs->get_fe(),
			   update_flags);

	::dealii::DataOutBase::Patch<DH::dimension, DH::space_dimension> sample_patch;
	sample_patch.n_subdivisions = n_subdivisions;
	sample_patch.data.reinit (n_datasets, patch_points.size());



	// now build the patches in parallel
	if (all_cells.size() > 0)
	dealii::WorkStream::run (&all_cells[0],
					 &all_cells[0]+all_cells.size(),
					 std_cxx1x::bind(&TEDataOut<dim,DH>::build_one_patch,
									 *this, std_cxx1x::_1, std_cxx1x::_2, std_cxx1x::_3,
									 curved_cell_region),
					 std_cxx1x::bind(&FEASolverNS::append_patch_to_list<dim,DH::space_dimension>,
									 std_cxx1x::_1, std_cxx1x::ref(this->patches)),
					 thread_data,
					 sample_patch);
}


template <int dim, class DH>
void
TEDataOut<dim,DH>::
build_one_patch (const std::pair<cell_iterator, unsigned int> *cell_and_index,
                 internal::DataOut::ParallelData<DH::dimension, DH::space_dimension> &data,
                 DataOutBase::Patch<DH::dimension, DH::space_dimension> &patch,
                 const CurvedCellRegion curved_cell_region)
{
  // use ucd_to_deal map as patch vertices
  // are in the old, unnatural ordering. if
  // the mapping does not preserve locations
  // (e.g. MappingQEulerian), we need to
  // compute the offset of the vertex for the
  // graphical output. Otherwise, we can just
  // use the vertex info.
  for (unsigned int vertex=0; vertex<GeometryInfo<DH::dimension>::vertices_per_cell; ++vertex)
    if (data.mapping_collection[0].preserves_vertex_locations())
      patch.vertices[vertex] = cell_and_index->first->vertex(vertex);
    else
      patch.vertices[vertex] = data.mapping_collection[0].transform_unit_to_real_cell
                               (cell_and_index->first,
                                GeometryInfo<DH::dimension>::unit_cell_vertex (vertex));

  if (data.n_datasets > 0)
    {
      data.x_fe_values.reinit (cell_and_index->first);
      const FEValues<DH::dimension,DH::space_dimension> &fe_patch_values
        = data.x_fe_values.get_present_fe_values ();

      const unsigned int n_q_points = fe_patch_values.n_quadrature_points;

      // depending on the requested output
      // of curved cells, if necessary
      // append the quadrature points to
      // the last rows of the patch.data
      // member. This is the case if we
      // want to produce curved cells at
      // the boundary and this cell
      // actually is at the boundary, or
      // else if we want to produce curved
      // cells everywhere
      //
      // note: a cell is *always* at
      // the boundary if dim<spacedim
      if (curved_cell_region==DataOut<dim, DH>::curved_inner_cells
          ||
          (curved_cell_region==DataOut<dim, DH>::curved_boundary
           &&
           (cell_and_index->first->at_boundary()
            ||
            (DH::dimension != DH::space_dimension))))
        {
          Assert(patch.space_dim==DH::space_dimension, ExcInternalError());
          const std::vector<Point<DH::space_dimension> > &q_points=fe_patch_values.get_quadrature_points();
          // resize the patch.data member
          // in order to have enough memory
          // for the quadrature points as
          // well
          patch.data.reinit (data.n_datasets+DH::space_dimension, n_q_points);
          // set the flag indicating that
          // for this cell the points are
          // explicitly given
          patch.points_are_available=true;
          // copy points to patch.data
          for (unsigned int i=0; i<DH::space_dimension; ++i)
            for (unsigned int q=0; q<n_q_points; ++q)
              patch.data(patch.data.size(0)-DH::space_dimension+i,q)=q_points[q][i];
        }
      else
        {
          patch.data.reinit(data.n_datasets, n_q_points);
          patch.points_are_available = false;
        }


      // counter for data records
      unsigned int offset=0;

      // first fill dof_data
      for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
        {
          const DataPostprocessor<DH::space_dimension> *postprocessor=this->dof_data[dataset]->postprocessor;
          Assert(postprocessor != 0, ExcMessage("TEDataOut requires the postprocessing class TEPostProcessor."));

          // Dynamically cast the postprocessor to the TMPostProcessor type, this allows the comput_derived_quantiites
          // to pass the material_id to the postprocessor
          const TEPostProcessor<DH::space_dimension> *tm_postprocessor = dynamic_cast<const TEPostProcessor<DH::space_dimension> * >(postprocessor);
          Assert(tm_postprocessor != 0, ExcMessage("Dynamic cast to TEPostProcessor, failed.  Be sure to use TEPostProcessor with TEDataOut"));

		  // Get material id from the cell
		  int mat_id = fe_patch_values.get_cell()->material_id();

		  // we have to postprocess the
		  // data, so determine, which
		  // fields have to be updated
		  const UpdateFlags update_flags=postprocessor->get_needed_update_flags();
		  if (data.n_components == 1)
			{
			  Assert(false, ExcMessage("TEDataOut should not be used for a single scalar data output. The point is to give stress output."));
			}
		  else
			{
			  // at each point there is
			  // a vector valued
			  // function and its
			  // derivative...
			  if (update_flags & update_values)
				this->dof_data[dataset]->get_function_values (fe_patch_values,
															  data.patch_values_system);
			  if (update_flags & update_gradients)
				this->dof_data[dataset]->get_function_gradients (fe_patch_values,
																 data.patch_gradients_system);
			  if (update_flags & update_hessians)
				this->dof_data[dataset]->get_function_hessians (fe_patch_values,
																data.patch_hessians_system);

			  if (update_flags & update_quadrature_points)
				data.patch_evaluation_points = fe_patch_values.get_quadrature_points();

			  std::vector<Point<DH::space_dimension> > dummy_normals;

			  tm_postprocessor->
			  compute_derived_quantities_vector(data.patch_values_system,
												data.patch_gradients_system,
												data.patch_hessians_system,
												dummy_normals,
												data.patch_evaluation_points,
												mat_id,
												data.postprocessed_values[dataset]);
			}

		  for (unsigned int q=0; q<n_q_points; ++q)
			for (unsigned int component=0;
				 component<this->dof_data[dataset]->n_output_variables;
				 ++component)
			  patch.data(offset+component,q)
				= data.postprocessed_values[dataset][q](component);
        }

      // then do the cell data. only
      // compute the number of a cell if
      // needed; also make sure that we
      // only access cell data if the
      // first_cell/next_cell functions
      // only return active cells
      if (this->cell_data.size() != 0)
        {
          Assert (!cell_and_index->first->has_children(), ExcNotImplemented());

          for (unsigned int dataset=0; dataset<this->cell_data.size(); ++dataset)
            {
              const double value
                = this->cell_data[dataset]->get_cell_data_value (cell_and_index->second);
              for (unsigned int q=0; q<n_q_points; ++q)
                patch.data(offset+dataset,q) = value;
            }
        }
    }


  for (unsigned int f=0; f<GeometryInfo<DH::dimension>::faces_per_cell; ++f)
    {
      // let's look up whether
      // the neighbor behind that
      // face is noted in the
      // table of cells which we
      // treat. this can only
      // happen if the neighbor
      // exists, and is on the
      // same level as this cell,
      // but it may also happen
      // that the neighbor is not
      // a member of the range of
      // cells over which we
      // loop, in which case the
      // respective entry in the
      // cell_to_patch_index_map
      // will have the value
      // no_neighbor. (note that
      // since we allocated only
      // as much space in this
      // array as the maximum
      // index of the cells we
      // loop over, not every
      // neighbor may have its
      // space in it, so we have
      // to assume that it is
      // extended by values
      // no_neighbor)
      if (cell_and_index->first->at_boundary(f)
          ||
          (cell_and_index->first->neighbor(f)->level() != cell_and_index->first->level()))
        {
          patch.neighbors[f] = numbers::invalid_unsigned_int;
          continue;
        }

      const cell_iterator neighbor = cell_and_index->first->neighbor(f);
      Assert (static_cast<unsigned int>(neighbor->level()) <
              data.cell_to_patch_index_map->size(),
              ExcInternalError());
      if ((static_cast<unsigned int>(neighbor->index()) >=
           (*data.cell_to_patch_index_map)[neighbor->level()].size())
          ||
          ((*data.cell_to_patch_index_map)[neighbor->level()][neighbor->index()]
           ==
           dealii::DataOutBase::Patch<DH::dimension>::no_neighbor))
        {
          patch.neighbors[f] = numbers::invalid_unsigned_int;
          continue;
        }

      // now, there is a
      // neighbor, so get its
      // patch number and set it
      // for the neighbor index
      patch.neighbors[f]
        = (*data.cell_to_patch_index_map)[neighbor->level()][neighbor->index()];
    }
}

template <int dim, class DH>
typename DataOut<dim,DH>::cell_iterator
TEDataOut<dim,DH>::first_locally_owned_cell ()
{
  typename DataOut<dim,DH>::cell_iterator
  cell = this->dofs->begin_active ();

  // skip cells if the current one
  // has no children (is active) and
  // is a ghost or artificial cell
  while ((cell != this->dofs->end()) &&
         (cell->has_children() == false) &&
         !cell->is_locally_owned())
    cell = DataOut<dim,DH>::next_cell(cell);

  return cell;
}



template <int dim, class DH>
typename DataOut<dim,DH>::cell_iterator
TEDataOut<dim,DH>::next_locally_owned_cell (const typename DataOut<dim,DH>::cell_iterator &old_cell)
{
  typename DataOut<dim,DH>::cell_iterator
  cell = DataOut<dim,DH>::next_cell(old_cell);
  while ((cell != this->dofs->end()) &&
         (cell->has_children() == false) &&
         !cell->is_locally_owned())
    cell = DataOut<dim,DH>::next_cell(cell);
  return cell;
}

/**
 * In a WorkStream context, use
 * this function to append the
 * patch computed by the parallel
 * stage to the array of patches.
 */
template <int dim, int spacedim>
void
append_patch_to_list (const DataOutBase::Patch<dim,spacedim> &patch,
					  std::vector<DataOutBase::Patch<dim,spacedim> > &patches)
{
	patches.push_back (patch);
	patches.back().patch_index = patches.size()-1;
}

}

#endif /* TEDATAOUT_H_ */
