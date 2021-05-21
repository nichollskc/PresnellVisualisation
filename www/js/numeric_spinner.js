$( function() {
    $( ".spinner" ).spinner( {
        min: 1,
        max: 60,
    });

    // Make sure the spinner updates on arrow clicking too
   $('.ui-spinner-button').click(function() {
       $(this).siblings('input').change();
   });
});


var spinnerBinding = new Shiny.InputBinding();
$.extend(spinnerBinding, {
  find: function(scope) {
    return $(scope).find(".spinner");
  },
  getValue: function(el) {
    return parseInt($(el).spinner("value"));
  },
  setValue: function(el, value) {
    $(el).spinner("value");
  },
  subscribe: function(el, callback) {
    $(el).on("change.spinnerBinding", function(e) {
      callback();
    });
  },
  unsubscribe: function(el) {
    $(el).off(".spinnerBinding");
  }
});

Shiny.inputBindings.register(spinnerBinding);
