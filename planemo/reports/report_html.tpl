<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{{ title }}</title>

    <!-- Bootstrap -->
    {{ bootstrap_style }}
    {{ custom_style }}

    <!-- HTML5 shim and Respond.js for IE8 support of HTML5 elements and media queries -->
    <!-- WARNING: Respond.js doesn't work if you view the page via file:// -->
    <!--[if lt IE 9]>
      <script src="https://oss.maxcdn.com/html5shiv/3.7.2/html5shiv.min.js"></script>
      <script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
    <![endif]-->

  </head>
  <body>

    <nav class="navbar navbar-inverse navbar-fixed-top" role="navigation">
      <div class="container-fluid">
        <div class="navbar-header">

          <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
            <span class="sr-only">Toggle navigation</span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="#">{{ title }}</a>
        </div>
        <div id="navbar" class="navbar-collapse collapse">
          <ul class="nav navbar-nav navbar-right">
            <li><a href="https://galaxyproject.org">Galaxy</a></li>
            <li><a href="https://planemo.readthedocs.org">Planemo</a></li>
          </ul>
          <div class="navbar-form navbar-right">
          </div>
        </div>
      </div>
    </nav>

    <div class="container-fluid">
      <div class="row">
        <div class="col-sm-3 col-md-2 sidebar">
          <ul class="nav nav-sidebar">
            <li><a href="#overview" class="text-success"><strong>Overview</strong></a></li>
          </ul>
          <ul class="nav nav-sidebar" id="nav-sidebar-tests">
          </ul>
        </div>
        <div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main">
          <!-- <h1 class="page-header">Tests</h1> -->
          <h2 id="overview">Overview</h2>
          <div id="overview-content"></div>
          <div class="progress">
          </div>
          <h2 id="tests">Tests</h2>
          <p>The remainder of this contains a description for each test executed to run these jobs.</p>
        </div>
      </div>
    </div>

    {{ jquery_script }}
    {{ bootstrap_script }}
    {{ custom_script }}
    <script>
      var testDataUrl = getUrlParameter("test_data_url");
      if(testDataUrl) {
      $.ajax(
        {'url': testDataUrl,
         'type': 'GET',
        }
        )
        .success(function(content) { renderTestResults( $.parseJSON(content) ); })
        .failure(function() { alert("Failed to load test data.")} );
      } else {
        var test_data = {{ json.dumps(raw_data) }};
        renderTestResults(test_data);
      }
    </script>
  </body>
</html>
